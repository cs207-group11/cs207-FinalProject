import parser

def test_ReactionParser_IOError():
    try:
        parser = ReactionParser("no_such_file")
    except IOError as err:
        assert(type(err) == IOError)

def test_ReactionParser_species():
    xml_filename = "rxns.xml"
    parser = ReactionParser(xml_filename)
    parser()
    assert parser.species == ['H', 'O', 'OH', 'H2', 'H2O', 'O2']
    
def test_ReactionParser_type():
    xml_filename = "rxns.xml"
    parser = ReactionParser(xml_filename)
    parser()
    assert parser.reaction_list[0].rxn_type == 'Elementary'
    
def test_ReactionParser_rate_coeffs_components():
    xml_filename = "rxns.xml"
    parser = ReactionParser(xml_filename)
    parser()
    assert parser.reaction_list[0].rate_coeffs_components == {'A': 35200000000.0, 'E': 71400.0}
    
def test_arr_A():
    try:
        xml_filename = "A_arr.xml"
        parser = ReactionParser(xml_filename)
        parser()
    except ValueError as err:
        assert(type(err) == ValueError)
        
def test_arr_E():
    try:
        xml_filename = "E_arr.xml"
        parser = ReactionParser(xml_filename)
        parser()
    except ValueError as err:
        assert(type(err) == ValueError)

def test_mod_arr_A():
    try:
        xml_filename = "A_mod_arr.xml"
        parser = ReactionParser(xml_filename)
        parser()
    except ValueError as err:
        assert(type(err) == ValueError)
        
def test_mod_arr_b():
    try:
        xml_filename = "b_mod_arr.xml"
        parser = ReactionParser(xml_filename)
        parser()
    except ValueError as err:
        assert(type(err) == ValueError)
        
def test_mod_arr_E():
    try:
        xml_filename = "E_mod_arr.xml"
        parser = ReactionParser(xml_filename)
        parser()
    except ValueError as err:
        assert(type(err) == ValueError)
        
def test_const_k():
    try:
        xml_filename = "k_const.xml"
        parser = ReactionParser(xml_filename)
        parser()
    except ValueError as err:
        assert(type(err) == ValueError)