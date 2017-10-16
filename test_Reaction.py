
"""Test module for Reaction.py"""

import pytest
from Reaction import *

def test_Reaction_special_fcns():
    test_reaction = Reaction(rxn_type="Elementary",
                             is_reversible=False,
                             rxn_equation="H2 + OH =] H2O + H",
                             rate_coeffs_components={'k': 0.5},
                             reactant_stoich_coeffs={'H2' :1, 'OH':1},
                            product_stoich_coeffs={'H2O' :1,  'H':1})
    info_reaction = "Reaction : H2 + OH =] H2O + H"
    assert str(test_reaction) == info_reaction
    assert len(test_reaction) == 4

def test_Reaction_set_bad_temperature():
    test_reaction = Reaction(rxn_type="Elementary",
                             is_reversible=False,
                             rxn_equation="H2 + OH =] H2O + H",
                             rate_coeffs_components={'k': 0.5},
                             reactant_stoich_coeffs={'H2' :1, 'OH':1},
                            product_stoich_coeffs={'H2O' :1,  'H':1})
    with pytest.raises(ValueError):
        test_reaction.set_temperature(0)


# def test_Reaction_set_bad_concentration():
#     test_reaction = Reaction(rxn_type="Elementary",
#                              is_reversible=False,
#                              rxn_equation="H2 + OH =] H2O + H",
#                              rate_coeffs_components={'k': 0.5},
#                              reactant_stoich_coeffs={'H2' :1, 'OH':1},
#                             product_stoich_coeffs={'H2O' :1,  'H':1})
#     with pytest.raises(ValueError):
#         pass

#     pass # TODO!
