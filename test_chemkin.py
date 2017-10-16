
"""Test module for chemkin."""

import numpy
import pytest
from chemkin import *

@pytest.fixture
def test_base_reaction():
    """Returns a valid reaction (from rxns.xml)"""
    return Reaction(rxn_type="Elementary",
                    is_reversible=False,
                    rxn_equation="H2 + OH =] H2O + H",
                    rate_coeffs_components={'k': 10},
                    reactant_stoich_coeffs={'H2' :1, 'OH':1},
                    product_stoich_coeffs={'H2O' :1, 'H':1})

def test_Reaction_special_fcns(test_base_reaction):
    """Test special functions for Reaction object"""
    info_reaction = "Reaction : H2 + OH =] H2O + H"
    assert str(test_base_reaction) == info_reaction
    assert len(test_base_reaction) == 4

def test_Reaction_get_unique_species(test_base_reaction):
    """Test get_unique_species functio (used to get
    values of species_list)."""
    expected = ['H2', 'H', 'H2O', 'OH']
    print test_base_reaction.species_list
    assert test_base_reaction.species_list == expected

def test_Reaction_set_zero_temperature(test_base_reaction):
    """Test setting reaction temperature to absolute 0..."""
    with pytest.raises(ValueError):
        test_base_reaction.set_temperature(0)

def test_Reaction_set_neg_temperature(test_base_reaction):
    """Test setting reaction temperature to a negative value"""
    with pytest.raises(ValueError):
        test_base_reaction.set_temperature(-100)

def test_Reaction_set_valid_temperature(test_base_reaction):
    """Test setting reaction temperature to absolute 0..."""
    test_base_reaction.set_temperature(10)
    assert test_base_reaction.temperature == 10

def test_Reaction_order_dictionaries(test_base_reaction):
    """Test ordering of reaction"""
    expected = [1, 4, 3, 2]
    test_list = test_base_reaction.order_dictionaries({'H2':1, 'OH':2, 'H2O':3, 'H':4})
    assert test_list == expected

def test_Reaction_set_concentrations(test_base_reaction):
    """Test setting reaction with valid concentrations"""
    expected = [1, 4, 3, 2]
    test_base_reaction.set_concentrations({'H2':1, 'OH':2, 'H2O':3, 'H':4})
    assert test_base_reaction.concentrations == expected

def test_Reaction_set_neg_concentrations(test_base_reaction):
    """Test setting reaction with negative concentrations"""
    with pytest.raises(ValueError):
        test_base_reaction.set_concentrations({'H2':1, 'OH':2, 'H2O':3, 'H':-4})

def test_Reaction_compute_reaction_rate_coeff(test_base_reaction):
    """Test calling compute_reaction_rate_coeff"""
    with pytest.raises(NotImplementedError):
        test_base_reaction.compute_reaction_rate_coeff()

def test_Reaction_compute_reaction_rate(test_base_reaction):
    """Test calling compute_reaction_rate"""
    with pytest.raises(NotImplementedError):
        test_base_reaction.compute_reaction_rate()