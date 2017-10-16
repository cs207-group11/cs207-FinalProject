import parser

def test_ReactionParser_species():
    xml_filename = "rxns.xml"
    parser = ReactionParser(xml_filename)
    parser()
    assert parser.species == ['H', 'O', 'OH', 'H2', 'O2']
    
def test_ReactionParser_type():
    xml_filename = "rxns.xml"
    parser = ReactionParser(xml_filename)
    parser()
    assert parser.reaction_list[0].rxn_type == 'Elementary'
    
def test_ReactionParser_rate_coeffs_components():
    xml_filename = "rxns.xml"
    parser = ReactionParser(xml_filename)
    parser()
    assert parser.reaction_list[0].rate_coeffs_components == {'A': 35200000000.0, 'E': 71400.0, 'b': -0.7}