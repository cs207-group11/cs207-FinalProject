from ReactionParser import * 

def test_ReactionParser_IOError():
    """Test when wrong xml is passed to parser"""
    try:
        parser = ReactionParser("no_such_file")
    except IOError as err:
        assert(type(err) == IOError)

def test_ReactionParser_species():
    """Test when reaction rate coefficient is modified
    Arrhenius but R is changed by user"""
    xml_filename = "rxns.xml"
    parser = ReactionParser(xml_filename)
    assert parser.species == {'H':None, 'O':None, 'OH':None, 'H2':None, 'H2O':None, 'O2':None}
    
def test_ReactionParser_type():
    """Test get_rxn_type() for an elementary reaction."""
    xml_filename = "rxns.xml"
    parser = ReactionParser(xml_filename)
    assert parser.reaction_list[0].rxn_type == 'Elementary'
    
def test_ReactionParser_rate_coeffs_components():
    """Test get_rate_coeffs_components for reaction 1."""
    xml_filename = "rxns.xml"
    parser = ReactionParser(xml_filename)
    assert parser.reaction_list[0].rate_coeffs_components == {'A': 35200000000.0, 'E': 71400.0}
    
def test_ReactionParser_is_reversible():
    """Test get_is_reversible for reaction irreversible reaction."""
    xml_filename = "rxns.xml"
    parser = ReactionParser(xml_filename)
    assert parser.reaction_list[0].is_reversible == False

def test_ReactionParser_rxn_equation():
    """Test get_rxn_equation for reaction 1."""
    xml_filename = "rxns.xml"
    parser = ReactionParser(xml_filename)
    assert parser.reaction_list[0].rxn_equation == 'H + O2 =] OH + O'
    
def test_ReactionParser_reactant_stoich_coeffs():
    """Test get_reactant_stoich_coeffs for reaction 1."""
    xml_filename = "rxns.xml"
    parser = ReactionParser(xml_filename)
    assert parser.reaction_list[0].reactant_stoich_coeffs == {'H': 1, 'H2': 0, 'H2O': 0, 'O': 0, 'O2': 1, 'OH': 0}

def test_ReactionParser_product_stoich_coeffs():
    """Test get_product_stoich_coeffs for reaction 1."""
    xml_filename = "rxns.xml"
    parser = ReactionParser(xml_filename)
    assert parser.reaction_list[0].product_stoich_coeffs == {'H': 0, 'H2': 0, 'H2O': 0, 'O': 1, 'O2': 0, 'OH': 1}
  
def test_unrecognizable_rxn():
    '''TEMPORARY'''
    """Test parser for reversible reaction. before milestone 2"""
    try:
        xml_filename = "unrecognized_rxn.xml"
        parser = ReactionParser(xml_filename)
    except NotImplementedError as err:
        assert(type(err) == NotImplementedError)

def test_arr_A():
    try:
        xml_filename = "A_arr.xml"
        parser = ReactionParser(xml_filename)
    except ValueError as err:
        assert(type(err) == ValueError)
        
def test_arr_E():
    try:
        xml_filename = "E_arr.xml"
        parser = ReactionParser(xml_filename)
    except ValueError as err:
        assert(type(err) == ValueError)

def test_mod_arr_A():
    try:
        xml_filename = "A_mod_arr.xml"
        parser = ReactionParser(xml_filename)
    except ValueError as err:
        assert(type(err) == ValueError)
        
def test_mod_arr_b():
    try:
        xml_filename = "b_mod_arr.xml"
        parser = ReactionParser(xml_filename)
    except ValueError as err:
        assert(type(err) == ValueError)
        
def test_mod_arr_E():
    try:
        xml_filename = "E_mod_arr.xml"
        parser = ReactionParser(xml_filename)
    except ValueError as err:
        assert(type(err) == ValueError)
        
def test_const_k():
    try:
        xml_filename = "k_const.xml"
        parser = ReactionParser(xml_filename)
    except ValueError as err:
        assert(type(err) == ValueError)
        
def test_get_NASA_poly_coef():
    parser = ReactionParser("rxns_reversible.xml")
    assert all(parser.NASA_poly_coefs[0]["high"]) == all([  2.50000000e+00,   0.00000000e+00,   0.00000000e+00,
         0.00000000e+00,   0.00000000e+00,   2.54716270e+04, -4.60117638e-01])