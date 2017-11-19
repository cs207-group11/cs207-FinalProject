
"""Crude interface/script for code review."""

from Parser import *

#User Input : Parse from input XML File
xml_filename1 = "rxns.xml"
xml_filename2 = "rxns_reversible.xml"
xml_filename3 = "rxns_mixed.xml"
parser = ReactionParser(xml_filename3)
#parser1 = ReactionParser(xml_filename1)
reaction_list = parser.get_reaction_list()
concentration2 = ({'H':1, 'H2':2, 'H2O':0, 'H2O2':0, 'HO2':0, 'O':0, "O2":2,"OH":0})
concentration1 = {'H':1, 'O2':2, 'OH':1, 'O':4, 'H2O':0, 'H2':1}
concentration3 = {'H':1, 'O2':2, 'OH':1, 'O':4, 'H2O':0, 'H2':1}
concentration_old = ({})
rxnsys = ReactionSystem(parser.reaction_list, parser.NASA_poly_coefs, 100, concentration3)
#rxnsys2 = ReactionSystem(parser1.reaction_list, parser.NASA_poly_coefs, 100, concentration1)


#Sanity Check for Testing
# for rxn1 in rxnsys.reaction_list:
#     print (rxn1.compute_reaction_rate_coeff())
#     print (rxn1.compute_progress_rate())
#     print (rxn1.compute_reaction_rate())

#User Call : Compute Reaction Rate Coefficient (K) --> Next : Progress Rate (and/or) Reaction Rate.
rxnrate = rxnsys.get_reaction_rate()
#rxnrate2 = rxnsys2.get_reaction_rate()

print("Reaction rate: {}".format(rxnrate))
#print("Reaction rate: {}".format(rxnrate2))