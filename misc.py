from libsbml import *

def dec_to_bin_tup(dim, num):
    """Convert decimal to binary tuple

    Args:
        dim: Size of tuple.
        num: Integer.

    Returns:
        tuple(base_2): m-tuple representation of num in base-2

    """

    pos_val = num
    base_2 = []
    base_2.append(pos_val % 2)
    for _ in range(dim - 1):
        pos_val = pos_val // 2
        base_2.append(pos_val % 2)
    return tuple(base_2)

def sbml_to_dict(file_name):
    reader = SBMLReader()
    document = reader.readSBML(file_name)

    if not document.getNumErrors():
        model = document.getModel()

        species = dict()
        iter_dict = {'synthesis': [2*i for i in range(100)], 'degradation': [2*i+1 for i in range(100)]}

        # go over all genes except for '_void_'
        for i in range(model.getNumSpecies() - 1):
            species_id = model.getSpecies(i).getId()
            species[species_id] = {}

            # add synthesis reaction of gene
            species[species_id]['synthesis'] = {}
            species[species_id]['synthesis']['parameters'] = {}
            species[species_id]['synthesis']['modifiers'] = []

            # add synthesis reaction of gene
            species[species_id]['degradation'] = {}
            species[species_id]['degradation']['parameters'] = {}
            species[species_id]['degradation']['modifiers'] = []


            for reaction_type in ['synthesis', 'degradation']:
                modifiers_list = model.getReaction(iter_dict[reaction_type][i]).getListOfModifiers()
                for j in range(len(modifiers_list)):
                    species[species_id][reaction_type]['modifiers'].append(modifiers_list[j].getSpecies())

                param_list = model.getReaction(iter_dict[reaction_type][i]).getKineticLaw().getListOfParameters()
                for j in range(len(param_list)):
                    species[species_id][reaction_type]['parameters'][param_list[j].getId()] = param_list[j].getValue()
        return species
    else:
        print('smbl file could not be parsed')
        return dict()