from misc import *


class KineticGRN:
    
    def __init__(self, grn_config=None):
        if not grn_config:
            self.genes = []
        else:
            self.genes = self.init_genes(grn_config)
            
    def init_genes(self, genes_dict):
        """
        define Genes and RegulatoryModules based on the GRN configuration dict.
        """
        
        genes = list(genes_dict.keys())
        genes_list = []

        for gene in genes_dict:
            modules = []
            # make a copy of modifiers list so we don't mutate the original one
            modifiers = [i for i in genes_dict[gene]['synthesis']['modifiers']]
            # little trick for counting number of occurances of module components
            params = list(genes_dict[gene]['synthesis']['parameters'].keys())
            num_modules = sum('bindsAsComplex' in s.split('_')[0] for s in params)
            num_modifiers = len(modifiers)
            relative_activations = []

            # create module object for each module entry in the dictionary
            p = 1 # placeholder
            for i in range(1,num_modules+1):
                binds_as_complex = genes_dict[gene]['synthesis']['parameters']['bindsAsComplex_{}'.format(i)]
                num_activators = int(genes_dict[gene]['synthesis']['parameters']['numActivators_{}'.format(i)])
                num_deactivators = int(genes_dict[gene]['synthesis']['parameters']['numDeactivators_{}'.format(i)])

                # fetch relevant parameters
                k = []
                n = []
                for j in range(p, p + num_activators + num_deactivators):
                    k.append(genes_dict[gene]['synthesis']['parameters']['k_{}'.format(j)])
                    n.append(genes_dict[gene]['synthesis']['parameters']['n_{}'.format(j)])

                p += num_activators + num_deactivators

                # populate the module object
                m = KineticRegulatoryModule()
                m.binds_as_complex = bool(binds_as_complex)
                m.modifiers = [genes.index(modifiers.pop(0)) for _ in range(num_activators + num_deactivators)]
                m.num_activators = num_activators
                m.num_deactivators = num_deactivators
                m.k = k
                m.n = n
                modules.append(m)

            # fetch relative activation coefficients
            for i in range(2**len(modules)):
                relative_activations.append(genes_dict[gene]['synthesis']['parameters']['a_{}'.format(i)])

            # creating gene module and populating it with all necessary information
            max_transcription = genes_dict[gene]['synthesis']['parameters']['max']
            max_translation = genes_dict[gene]['synthesis']['parameters']['maxTranslation']
            delta_mRNA_deg = genes_dict[gene]['degradation']['parameters']['delta']
            delta_protein_deg = genes_dict[gene]['synthesis']['parameters']['deltaProtein']
            g = KineticGene(
                max_transcription,
                max_translation,
                delta_mRNA_deg,
                delta_protein_deg
            )
            g.relative_activations = relative_activations
            g.modules = modules
            genes_list.append(g)

        return genes_list

class KineticGene:
    """
    Gene within a GRN, which is associated to activating and repressing genes
    """
    
    def __init__(self, m_transcr, m_transl, d_mrna, d_prot):
        """
        m_trancr: maximum transcription rate
        m_transl: maximum translation rate
        d_mrna: mRNA degradation rate
        d_prot: protein degradation rate
        """

        self.max_transcription = m_transcr
        self.max_translation = m_transl
        self.delta_mRNA_deg = d_mrna
        self.delta_protein_deg = d_prot
        # list of relative activations for all possible TF binding configurations
        self.relative_activations = []
        self.wild_type_activations = [] # backup
        # list of transcription factors (TFs), also known as regulatory modules
        self.modules = []
        
        
    def input_func(self, y_vec):
        
        activation = 0
        m = []
        for module in self.modules:
            m.append(module.module_activation(y_vec))
        
        for i in range(len(self.relative_activations)):
            
            if not len(self.modules):
                activation = self.relative_activations[0]
                break
                
            # the i-th configuration of the TF binding states
            config = dec_to_bin_tup(len(self.modules), i)
            
            # compute probability that promoter is in state relative i-th activation
            prod = 1
            for j, indicator in enumerate(config):
                if indicator:
                    prod *= m[j]
                else:
                    prod *= 1 - m[j]
                    
            # overall gene activation
            activation += self.relative_activations[i] * prod

        return activation

    def perturb_activation(self, perturb):
        
        # make backup of original (wildtype) activations
        self.wild_type_activations = [i for i in self.relative_activations]

        # ensure perturbation allows activation to be in interval [0,1]
        if self.relative_activations[0] + perturb > 1:
            perturb = 1 - self.relative_activations[0]
        elif self.relative_activations[0] + perturb < 0:
            perturb = 0 - self.relative_activations[0]

        for i in range(len(self.relative_activations)):
            self.relative_activations[i] += perturb
            # truncate to [0 1]
            if self.relative_activations[i] < 0:
                self.relative_activations[i] = 0
            elif self.relative_activations[i] > 1:
                self.relative_activations[i] = 1

    def restore_activation(self):
        self.relative_activations = [i for i in self.wild_type_activations]
    

class KineticRegulatoryModule:
    """
    A module can be seen as a transcription factor (TF), which may be a complex of multiple genes interacting at the regulatory site
    """
    
    def __init__(self):
        # indicator for whether modifiers form a complex at the binding site
        self.binds_as_complex = False
        # number of modifiers that are activating
        self.num_activators = 0
        # number of modifiers that are deactivating or repressing
        self.num_deactivators = 0
        # the dissociation constants
        self.k = []
        # the hill coefficients
        self.n = []
        # the modifiers are actually indices so they can be fetched from y_vec
        self.modifiers = []
        
    def module_activation(self, y_vec):
        """
        Computes the module activation according to Hill's kinetic equation

        y: vector
        gene: the key
        modifiers: the genes that interact within the module
        """
        
        x = []
        for idx, i in enumerate(self.modifiers):
            x.append((y_vec[i]/self.k[idx])**self.n[idx])
        
        numer = 1
        for i in range(self.num_activators):
            numer *= x[i]
        
        denom = 1
        if self.binds_as_complex:
            # denom is 1 + prod_activ + prod_deactiv because of complex formation
            denom += numer
            
            if self.num_deactivators > 0:
                prod = 1
                for i in range(self.num_activators, len(self.modifiers)):
                    prod *= x[i]
                denom += prod
        else:
            for i in range(len(self.modifiers)):
                denom *= (1 + x[i])
        
        return numer / denom