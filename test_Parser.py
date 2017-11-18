
import pytest

from Parser import ReactionParser

# ======================= TESTS FOR REACTIONPARSER OBJECT ====================== #

def test_ReactionParser_file_not_found():
    """Test when xml file is nonexistent"""
    with pytest.raises(IOError):
        parser = ReactionParser("no_such_file")

def test_ReactionParser_species():
    """Test when reaction rate coefficient is modified
    Arrhenius but R is changed by user"""
    xml_filename = "rxns.xml"
    parser = ReactionParser(xml_filename)
    assert parser.get_species() == ({'H': None,'O': None, 'OH': None,
                                    'H2': None, 'H2O': None, 'O2': None})
    
def test_ReactionParser_type():
    """Test get_rxn_type() for an elementary reaction."""
    xml_filename = "rxns.xml"
    parser = ReactionParser(xml_filename)
    parser()
    assert parser.reaction_list[0].rxn_type == 'Elementary'
    
def test_ReactionParser_rate_coeffs_components():
    """Test get_rate_coeffs_components for reaction 1."""
    xml_filename = "rxns.xml"
    parser = ReactionParser(xml_filename)
    parser()
    assert (parser.reaction_list[0].rate_coeffs_components ==
            {'A': 35200000000.0, 'E': 71400.0})
    
def test_ReactionParser_is_reversible():
    """Test get_is_reversible for reaction irreversible reaction."""
    xml_filename = "rxns.xml"
    parser = ReactionParser(xml_filename)
    parser()
    assert parser.reaction_list[0].is_reversible == False

def test_ReactionParser_rxn_equation():
    """Test get_rxn_equation for reaction 1."""
    xml_filename = "rxns.xml"
    parser = ReactionParser(xml_filename)
    parser()
    assert parser.reaction_list[0].rxn_equation == 'H + O2 =] OH + O'
    
def test_ReactionParser_reactant_stoich_coeffs():
    """Test get_reactant_stoich_coeffs for reaction 1."""
    xml_filename = "rxns.xml"
    parser = ReactionParser(xml_filename)
    parser()
    assert (parser.reaction_list[0].reactant_stoich_coeffs ==
            {'H': 1, 'H2': 0, 'H2O': 0, 'O': 0, 'O2': 1, 'OH': 0})

def test_ReactionParser_product_stoich_coeffs():
    """Test get_product_stoich_coeffs for reaction 1."""
    xml_filename = "rxns.xml"
    parser = ReactionParser(xml_filename)
    parser()
    assert (parser.reaction_list[0].product_stoich_coeffs ==
            {'H': 0, 'H2': 0, 'H2O': 0, 'O': 1, 'O2': 0, 'OH': 1})

def test_arr_A():
    """Test when parameter A (for computing
    Arrrhenius reaction rate coeff)
    is missing from xml file"""
    with pytest.raises(ValueError):
        xml_filename = "A_arr.xml"
        parser = ReactionParser(xml_filename)
        parser()
        
def test_arr_E():
    """Test when parameter E (for computing
    Arrrhenius reaction rate coeff)
    is missing from xml file"""
    with pytest.raises(ValueError):
        xml_filename = "E_arr.xml"
        parser = ReactionParser(xml_filename)
        parser()

def test_mod_arr_A():
    """Test when parameter A (for computing
    modified Arrrhenius reaction rate coeff)
    is missing from xml file"""
    with pytest.raises(ValueError):
        xml_filename = "A_mod_arr.xml"
        parser = ReactionParser(xml_filename)
        parser()

def test_mod_arr_b():
    """Test when parameter b (for computing
    modified Arrrhenius reaction rate coeff)
    is missing from xml file"""
    with pytest.raises(ValueError):
        xml_filename = "b_mod_arr.xml"
        parser = ReactionParser(xml_filename)
        parser()
        
def test_mod_arr_E():
    """Test when parameter E (for computing
    modified Arrrhenius reaction rate coeff)
    is missing from xml file"""
    with pytest.raises(ValueError):
        xml_filename = "E_mod_arr.xml"
        parser = ReactionParser(xml_filename)
        parser()
        
def test_const_k():
    """Test when k (for computing
    constant reaction rate coeff)
    is missing from xml file"""
    with pytest.raises(ValueError):
        xml_filename = "k_const.xml"
        parser = ReactionParser(xml_filename)
        parser()