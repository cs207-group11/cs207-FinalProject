
"""Test module for Parser."""

import os
import pytest

from chemkinlib.config import DATA_DIRECTORY
from chemkinlib.utils import Parser

# ======================= TESTS FOR REACTIONPARSER OBJECT ====================== #

def test_RxnParser_file_not_found():
    """Test when xml file is nonexistent"""
    with pytest.raises(IOError):
        parser = Parser.ReactionParser("no_such_file")

def test_RxnParser_species():
    """Test when reaction rate coefficient is modified
    Arrhenius but R is changed by user"""
    xml_filename = os.path.join(DATA_DIRECTORY, "rxns.xml")
    parser = Parser.ReactionParser(xml_filename)
    assert parser.get_species() == ({'H': None,'O': None, 'OH': None,
                                    'H2': None, 'H2O': None, 'O2': None})
    
def test_RxnParser_type():
    """Test get_rxn_type() for an elementary reaction."""
    xml_filename = os.path.join(DATA_DIRECTORY, "rxns.xml")
    parser = Parser.ReactionParser(xml_filename)
    assert parser.reaction_list[0].rxn_type == 'Elementary'
    
def test_RxnParser_rate_coeffs_components():
    """Test get_rate_coeffs_components for reaction 1."""
    xml_filename = os.path.join(DATA_DIRECTORY, "rxns.xml")
    parser = Parser.ReactionParser(xml_filename)
    assert (parser.reaction_list[0].rate_coeffs_components ==
            {'A': 35200000000.0, 'E': 71400.0})
    
def test_RxnParser_is_reversible():
    """Test get_is_reversible for reaction irreversible reaction."""
    xml_filename = os.path.join(DATA_DIRECTORY, "rxns.xml")
    parser = Parser.ReactionParser(xml_filename)
    assert parser.reaction_list[0].is_reversible == False

def test_RxnParser_rxn_equation():
    """Test get_rxn_equation for reaction 1."""
    xml_filename = os.path.join(DATA_DIRECTORY, "rxns.xml")
    parser = Parser.ReactionParser(xml_filename)
    assert parser.reaction_list[0].rxn_equation == 'H + O2 =] OH + O'
    
def test_RxnParser_reactant_stoich_coeffs():
    """Test get_reactant_stoich_coeffs for reaction 1."""
    xml_filename = os.path.join(DATA_DIRECTORY, "rxns.xml")
    parser = Parser.ReactionParser(xml_filename)
    assert (parser.reaction_list[0].reactant_stoich_coeffs ==
            {'H': 1, 'H2': 0, 'H2O': 0, 'O': 0, 'O2': 1, 'OH': 0})

def test_RxnParser_product_stoich_coeffs():
    """Test get_product_stoich_coeffs for reaction 1."""
    xml_filename = os.path.join(DATA_DIRECTORY, "rxns.xml")
    parser = Parser.ReactionParser(xml_filename)
    assert (parser.reaction_list[0].product_stoich_coeffs ==
            {'H': 0, 'H2': 0, 'H2O': 0, 'O': 1, 'O2': 0, 'OH': 1})

def test_arr_A():
    """Test when parameter A (for computing
    Arrrhenius reaction rate coeff)
    is missing from xml file"""
    with pytest.raises(ValueError):
        xml_filename = os.path.join(DATA_DIRECTORY, "A_arr.xml")
        parser = Parser.ReactionParser(xml_filename)
        
def test_arr_E():
    """Test when parameter E (for computing
    Arrrhenius reaction rate coeff)
    is missing from xml file"""
    with pytest.raises(ValueError):
        xml_filename = os.path.join(DATA_DIRECTORY, "E_arr.xml")
        parser = Parser.ReactionParser(xml_filename)

def test_mod_arr_A():
    """Test when parameter A (for computing
    modified Arrrhenius reaction rate coeff)
    is missing from xml file"""
    with pytest.raises(ValueError):
        xml_filename = os.path.join(DATA_DIRECTORY, "A_mod_arr.xml")
        parser = Parser.ReactionParser(xml_filename)

def test_mod_arr_b():
    """Test when parameter b (for computing
    modified Arrrhenius reaction rate coeff)
    is missing from xml file"""
    with pytest.raises(ValueError):
        xml_filename = os.path.join(DATA_DIRECTORY, "b_mod_arr.xml")
        parser = Parser.ReactionParser(xml_filename)
        
def test_mod_arr_E():
    """Test when parameter E (for computing
    modified Arrrhenius reaction rate coeff)
    is missing from xml file"""
    with pytest.raises(ValueError):
        xml_filename = os.path.join(DATA_DIRECTORY, "E_mod_arr.xml")
        parser = Parser.ReactionParser(xml_filename)
        
def test_const_k():
    """Test when k (for computing
    constant reaction rate coeff)
    is missing from xml file"""
    with pytest.raises(ValueError):
        xml_filename = os.path.join(DATA_DIRECTORY, "k_const.xml")
        parser = Parser.ReactionParser(xml_filename)

def test_get_coeffs_invalid_temp_range():
    """Test when invalid temperature range inputed"""
    with pytest.raises(ValueError):
        xml_filename = os.path.join(DATA_DIRECTORY, "rxns.xml")
        parser = Parser.ReactionParser(xml_filename)
        parser.get_coeffs('H', 'invalid_range')
