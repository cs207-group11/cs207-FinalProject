
# """Crude interface/script for code review."""

from Parser import *
from chemkin import *
import copy

xml_filename = "rxns_reversible.xml"
parser = ReactionParser(xml_filename)
concentration = ({'H':1, 'H2':1, 'H2O':1, 'H2O2':1, 'HO2':1, 'O':1, "O2":1, "OH":1})
temperature = 1000
rxnsys = ReactionSystem(parser.reaction_list,
                        parser.NASA_poly_coefs,
                        temperature,
                        concentration)



# xml_filename = "rxnset_long.xml"
# parser = ReactionParser(xml_filename)
# concentration = ({'H':1, 'H2':1, 'H2O':1, 'H2O2':1, 'HO2':1, 'O':1, "O2":1, "OH":1})
# temperature = 1000
# rxnsys = ReactionSystem(parser.reaction_list,
#                         parser.NASA_poly_coefs,
#                         temperature,
#                         concentration)




# #Sanity Check for Testing
# rates = []
# flag = 1
# for rxn in rxnsys.reaction_list:
#     #print(rxn, type(rxn))
#     print(rxn.species_list)
#     list_species_ordered = list(rxn.species_list)
#     rate = rxn.compute_reaction_rate()
#     rates.append(rate)

# rates = numpy.array(rates)
# print(rates.shape)

# rates_rxn_sys = numpy.sum(rates, axis=0)
# print(rates_rxn_sys)

#3print(rxnsys.involved_species)

rxnrates_dict = rxnsys.sort_reaction_rates()
# print(rxnrates_dict)

for k, v in rxnrates_dict.items():
    print("{0} \t {1:e}".format(k, v))
