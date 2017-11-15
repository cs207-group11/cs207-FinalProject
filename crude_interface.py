
"""Crude interface/script for code review."""

from chemkin import *

#User Input : Parse from input XML File
xml_filename = "rxns.xml"
parser = ReactionParser(xml_filename)
parser()
rxnsys = ReactionSystem(parser.reaction_list)


#User Input : Enter all Specie concentrations, order does not matter.
rxnsys.set_concentrations({'H':1, 'O2':2, 'OH':1, 'O':4, 'H2O':0, 'H2':1})

#User Input : Enter Temperature.
rxnsys.set_temperature(100)

#Sanity Check for Testing
for rxn1 in rxnsys.reaction_list:
    print (rxn1)
    print (rxn1.compute_reaction_rate_coeff())
    print (rxn1.compute_progress_rate())
    print (rxn1.compute_reaction_rate())

#User Call : Compute Reaction Rate Coefficient (K) --> Next : Progress Rate (and/or) Reaction Rate.
rxnrate = rxnsys.get_reaction_rate()

print("Reaction rate: {}".format(rxnrate))
