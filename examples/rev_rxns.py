
"""Example of reversible reaction."""

import os

from chemkinlib.utils import Parser
from chemkinlib.reactions import ReactionSystems
from chemkinlib.config import DATA_DIRECTORY

# USER INPUT: reaction (xml) file
xml_filename = os.path.join(DATA_DIRECTORY, "rxns_reversible.xml")

parser = Parser.ReactionParser(xml_filename)

# USER INPUTS (concentrations and temperatures)
concentration = ({'H':1, 'H2':1, 'H2O':1, 'H2O2':1, 'HO2':1, 'O':1, "O2":1, "OH":1})
temperature = 1000


# Set up reaction system
rxnsys = ReactionSystems.ReactionSystem(parser.reaction_list,
                                  parser.NASA_poly_coefs,
                                  temperature,
                                  concentration)

# Compute and sort reaction rates
rxnrates_dict = rxnsys.sort_reaction_rates()

# display reaction rates by species
for k, v in rxnrates_dict.items():
    print("d[{0}]/dt : \t {1:e}".format(k, v))
