
"""Crude interface/script for code review."""

from chemkin import *

### [USER INPUT REQUIRED] Parse from input xml file
xml_filename = "rxns.xml"
parser = ReactionParser(xml_filename)
parser()
rxn1 = parser.reaction_list[0]


print rxn1
# print rxn1.species_list
# print rxn1.reactant_stoich_coeffs
# print rxn1.product_stoich_coeffs
# print rxn1.rate_coeffs_components


### [USER INPUT REQUIRED] Be sure to enter in cocentrations for the right species!
rxn1.set_concentrations({'H':1, 'O2':2, 'OH':0, 'O':0, 'H2O':0, 'H2':0})

### [USER INPUT REQUIRED]
rxn1.set_temperature(100)


k = rxn1.compute_reaction_rate_coeff()
omega = rxn1.compute_progress_rate()
rxnrate = rxn1.compute_reaction_rate()

print("Reaction rate coefficient (k) : {}".format(k))
print("Reaction progress rate: {}".format(omega))
print("Reaction rate: {}".format(rxnrate))
