
"""Crude interface/script for code review."""

from chemkin import *

#User Input : Parse from input XML File
xml_filename = "rxns.xml"
parser = ReactionParser(xml_filename)
parser()
rxn1 = parser.reaction_list[0]

#Sanity Check for Testing
print(rxn1)
# print rxn1.species_list
# print rxn1.reactant_stoich_coeffs
# print rxn1.product_stoich_coeffs
# print rxn1.rate_coeffs_components


#User Input : Enter all Specie concentrations, order does not matter.
rxn1.set_concentrations({'H':1, 'O2':2, 'OH':0, 'O':0, 'H2O':0, 'H2':0})

#User Input : Enter Temperature.
rxn1.set_temperature(100)

#User Call : Compute Reaction Rate Coefficient (K) --> Next : Progress Rate (and/or) Reaction Rate.
k = rxn1.compute_reaction_rate_coeff()
omega = rxn1.compute_progress_rate()
rxnrate = rxn1.compute_reaction_rate()

print("Reaction rate coefficient (k) : {}".format(k))
print("Reaction progress rate: {}".format(omega))
print("Reaction rate: {}".format(rxnrate))
