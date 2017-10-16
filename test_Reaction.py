
"""Test module for Reaction.py"""

import pytest
from Reaction import *

@pytest.fixture
def test_reaction_1():
    """Returns a valid reaction (from rxns.xml)"""
    return Reaction(rxn_type="Elementary",
                    is_reversible=False,
                    rxn_equation="H2 + OH =] H2O + H",
                    rate_coeffs_components={'k': 0.5},
                    reactant_stoich_coeffs={'H2' :1, 'OH':1},
                    product_stoich_coeffs={'H2O' :1,  'H':1})

def test_Reaction_special_fcns(test_reaction_1):
    """Test special functions for Reaction object"""
    info_reaction = "Reaction : H2 + OH =] H2O + H"
    assert str(test_reaction_1) == info_reaction
    assert len(test_reaction_1) == 4

def test_Reaction_set_0_temperature(test_reaction_1):
    """Test setting reaction temperature to absolute 0..."""
    with pytest.raises(ValueError):
        test_reaction_1.set_temperature(0)

def test_Reaction_set_neg_temperature(test_reaction_1):
    """Test setting reaction temperature to a negative value"""
    with pytest.raises(ValueError):
        test_reaction_1.set_temperature(-100)

# def test_Reaction_set_bad_concentration(test_reaction_2):
#     with pytest.raises(ValueError):
#         pass

#     pass # TODO!
