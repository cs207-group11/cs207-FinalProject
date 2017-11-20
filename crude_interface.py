
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



rxnrates_dict = rxnsys.sort_reaction_rates()

for k, v in rxnrates_dict.items():
    print("d[{0}]/dt : \t {1:e}".format(k, v))
