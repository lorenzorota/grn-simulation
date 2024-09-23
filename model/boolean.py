from misc import *


class BooleanGRN:
    """
    Stores additional information about the Boolean network, such as modules, modifiers, relative activations
    """
    def __init__(self, grn_config=None):
        if not grn_config:
            self.genes = []
        else:
            self.genes = self.init_genes(grn_config)
            
    def init_genes(self, genes_dict):
        """
        define Genes and RegulatoryModules based on the GRN configuration dict.
        """
        
        genes = dict()

        for gene in genes_dict:
            modules = []
            # make a copy of modifiers list so we don't mutate the original one
            modifiers = [i for i in genes_dict[gene]['synthesis']['modifiers']]
            # little trick for counting number of occurances of module components
            params = list(genes_dict[gene]['synthesis']['parameters'].keys())
            num_modules = sum('bindsAsComplex' in s.split('_')[0] for s in params)
            num_modifiers = len(modifiers)
            relative_activations = []
            boolean_activations = []
            
            # create module object for each module entry in the dictionary
            for i in range(1,num_modules+1):
                binds_as_complex = genes_dict[gene]['synthesis']['parameters']['bindsAsComplex_{}'.format(i)]
                num_activators = int(genes_dict[gene]['synthesis']['parameters']['numActivators_{}'.format(i)])
                num_deactivators = int(genes_dict[gene]['synthesis']['parameters']['numDeactivators_{}'.format(i)])

                # populate the module object
                m = BooleanRegulatoryModule()
                m.binds_as_complex = bool(binds_as_complex)
                m.modifiers = [modifiers.pop(0) for _ in range(num_activators + num_deactivators)]
                m.num_activators = num_activators
                m.num_deactivators = num_deactivators
                modules.append(m)

            # fetch relative activation coefficients and discretize them
            for i in range(2**len(modules)):
                activation = genes_dict[gene]['synthesis']['parameters']['a_{}'.format(i)]
                bool_activation = True if activation >= 0.5 else False
                relative_activations.append(activation)
                boolean_activations.append(bool_activation)

            # creating gene and populating it with all necessary information
            g = BooleanGene(gene)
            g.relative_activations = relative_activations
            g.boolean_activations = boolean_activations
            g.modules = modules
            genes[gene] = g

        return genes

class BooleanGene:
    """
    Gene within a Boolean GRN, which is associated to activating and repressing genes
    """
    
    def __init__(self, gene):
        """
        """
        
        self.gene = gene
        # list of relative activations for all possible TF binding configurations
        self.relative_activations = []
        self.wild_type_activations = [] # backup
        self.boolean_activations = [] # after discretizing
        # list of transcription factors (TFs), also known as regulatory modules
        self.modules = []

    def input_function(self, G):
        """ Computes the relative activation of nodes """
        
        # initialize list of modules
        m = []
        
        # first compute module binding state
        for module in self.modules:
            if module.binds_as_complex:
                activators = True
                repressors = True
                
                # first compute product of activators
                for i in range(module.num_activators):
                    modif = module.modifiers[i]
                    action = G[modif][self.gene]['weight']
                    state = G.nodes[modif]['state']
                    repressors = repressors & (action & state)
                    
                # next compute product of deactivators
                for i in range(module.num_activators, len(module.modifiers)):
                    modif = module.modifiers[i]
                    action = G[modif][self.gene]['weight']
                    state = G.nodes[modif]['state']
                    activators = activators & (action & state)
                
                # because we have True xor (action and state) xor ...
                sum_ = True ^ activators ^ repressors
                m.append(sum_)
            else:
                prod_ = True
                for modif in module.modifiers:
                    action = G[modif][self.gene]['weight']
                    state = G.nodes[modif]['state']
                    prod_ = prod_ & (True ^ (action & state)) # (True xor (action and state)) and ...
                m.append(prod_)
        
        # compute final decision based on relative activations of module binding state configurations
        out = False # neutral value in xor chain
        for i in range(len(self.boolean_activations)):

            # activation has no effect if node is not regulated
            if not len(self.modules):
                return G.nodes[self.gene]['state']

            # the i-th configuration of the TF binding states
            config = dec_to_bin_tup(len(self.modules), i)

            prod_ = True
            for j, indicator in enumerate(config):
                if indicator:
                    prod_ = prod_ & m[j]
                else:
                    prod_ = prod_ & (not m[j])

            # overall gene activation
            out = out ^ (self.boolean_activations[i] & prod_)
        return out

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
        
        self.boolean_activations =  [True if i >= 0.5 else False for i in self.relative_activations]

    def restore_activation(self):
        self.relative_activations = [i for i in self.wild_type_activations]
        self.boolean_activations =  [True if i >= 0.5 else False for i in self.relative_activations]

        
class BooleanRegulatoryModule:
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
        # the modifiers are actually indices so they can be fetched from y_vec
        self.modifiers = []