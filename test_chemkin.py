"""Test module for chemkin."""

import numpy
import pytest
import warnings
from chemkin import *
from Parser import ReactionParser

# Treat warnings like errors (for testing purposes)
warnings.simplefilter("error")


# ======================= TESTS FOR REACTIONSYSTEM OBJECT ====================== #

@pytest.fixture
def test_rxn_sys():
    """Returns a valid reaction system"""
    xml_filename = "rxn.xml"
    parser = ReactionParser(xml_filename)
    temp = 500 # "low" temperature range in NASA coeffs database
    concentrations = {'H':1, 'O2':1, 'H2O':1}
    rxnsys = ReactionSystem(parser.reaction_list, parser.NASA_poly_coefs, temp, concentrations)
    return rxnsys

def test_rxn_sys_invalid_temperature():
    """Tests setting up reaction system with invalid temperatures."""
    xml_filename = "rxns_mixed.xml"
    parser = ReactionParser(xml_filename)
    concentrations = {'H':1, 'O2':2, 'OH':1, 'O':4, 'H2O':0, 'H2':1}
    temp = 0
    with pytest.raises(ValueError):
        rxnsys = ReactionSystem(parser.reaction_list, parser.NASA_poly_coefs, temp, concentrations)
    temp = -100
    with pytest.raises(ValueError):
        rxnsys = ReactionSystem(parser.reaction_list, parser.NASA_poly_coefs, temp, concentrations)


def test_rxn_sys_get_reaction_rate_for_1_rxn(test_rxn_sys):
    """Tests function to get reaction rate for a given system of reactions (just 1 reaction)."""
    assert test_rxn_sys.get_reaction_rate() == -15.

def test_rxn_sys_get_reaction_rate_for_3_rxns():
    """Tests function to get reaction rate for a given system of reactions (more than 1 reaction)."""
    xml_filename = "rxnsys.xml"
    parser = ReactionParser(xml_filename)
    temp = 10
    concentrations = {'H':1, 'O2':1, 'OH':1, 'O':1, 'H2O':1, 'H2':1}
    rxnsys = ReactionSystem(parser.reaction_list, parser.NASA_poly_coefs, temp, concentrations)
    assert rxnsys.get_reaction_rate() == -10.

def test_rxn_sys_get_lowT_nasa_matrix(test_rxn_sys):
    """Tests function to fetch NASA coefficients of appropriate T and appropriate species in reaction."""
    # Order of low T range NASA coefficients: H, H2O, O2
    expected_nasa = {'H': numpy.array([2.50000000e+00, 0.00000000e+00,
                                       0.00000000e+00, 0.00000000e+00, 
                                       0.00000000e+00, 2.54716270e+04, -4.60117608e-01]),
                    'H2O': numpy.array([3.38684249e+00, 3.47498246e-03,
                                        -6.35469633e-06, 6.96858127e-09,
                                        -2.50658847e-12, -3.02081133e+04, 2.59023285e+00]),
                    'O2': numpy.array([3.21293640e+00, 1.12748635e-03,
                                       -5.75615047e-07, 1.31387723e-09,
                                       -8.76855392e-13, -1.00524902e+03, 6.03473759e+00])}

    
    # expected_nasa = numpy.array([[2.50000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 
    #                             0.00000000e+00, 2.54716270e+04, -4.60117608e-01],
    #                             [3.38684249e+00, 3.47498246e-03, -6.35469633e-06, 6.96858127e-09,
    #                             -2.50658847e-12, -3.02081133e+04, 2.59023285e+00],
    #                             [3.21293640e+00, 1.12748635e-03, -5.75615047e-07, 1.31387723e-09,
    #                             -8.76855392e-13, -1.00524902e+03, 6.03473759e+00]])
    assert (numpy.isclose(test_rxn_sys.NASA_matrix['H2O'], expected_nasa['H2O'])).all()
    assert (numpy.isclose(test_rxn_sys.NASA_matrix['O2'], expected_nasa['O2'])).all()
    assert (numpy.isclose(test_rxn_sys.NASA_matrix['H'], expected_nasa['H'])).all()


def test_rxn_sys_get_highT_nasa_matrix():
    """Tests function to fetch NASA coefficients of appropriate T and appropriate species in reaction."""
    xml_filename = "rxn.xml"
    parser = ReactionParser(xml_filename)
    temp = 5000 # "high" temperature range in NASA coeffs database
    concentrations = {'H':1, 'O2':1, 'H2O':1}
    rxnsys = ReactionSystem(parser.reaction_list, parser.NASA_poly_coefs, temp, concentrations)
    # Order of high T range NASA coefficients: H, H2O, O2
    expected_nasa = {'O2': numpy.array([3.69757819e+00, 6.13519689e-04,
                                       -1.25884199e-07, 1.77528148e-11,
                                       -1.13643531e-15,  -1.23393018e+03, 3.18916559e+00]),
                    'H2O': numpy.array([2.67214561e+00, 3.05629289e-03,
                                       -8.73026011e-07, 1.20099639e-10,
                                       -6.39161787e-15,  -2.98992090e+04, 6.86281681e+00]),
                    'H': numpy.array([2.50000000e+00, 0.00000000e+00,
                                     0.00000000e+00, 0.00000000e+00,
                                     0.00000000e+00, 2.54716270e+04, -4.60117638e-01])}
    # expected_nasa = numpy.array([[2.50000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 
    #                             0.00000000e+00, 2.54716270e+04, -4.60117608e-01],
    #                             [2.67214561e+00, 3.05629289e-03, -8.73026011e-07, 1.20099639e-10,
    #                             -6.39161787e-15, -2.98992090e+04, 6.86281681e+00],
    #                             [3.69757819e+00, 6.13519689e-04, -1.25884199e-07, 1.77528148e-11,
    #                             -1.13643531e-15, -1.23393018e+03, 3.18916559e+00]])
    #assert numpy.isclose(rxnsys.NASA_matrix, expected_nasa).all()
    assert (numpy.isclose(rxnsys.NASA_matrix['H2O'], expected_nasa['H2O'])).all()
    assert (numpy.isclose(rxnsys.NASA_matrix['O2'], expected_nasa['O2'])).all()
    assert (numpy.isclose(rxnsys.NASA_matrix['H'], expected_nasa['H'])).all()


def test_rxn_sys_rev_reaction():
    """Tests setting up reaction system with reversible reaction."""
    xml_filename = "rev_rxn.xml"
    parser = ReactionParser(xml_filename)
    temp = 500 # "low" temperature range in NASA coeffs database
    concentrations = {'H':1, 'O2':1, 'H2O':1}
    rxnsys = ReactionSystem(parser.reaction_list, parser.NASA_poly_coefs, temp, concentrations)
    # expected_nasa = numpy.array([[2.50000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 
    #                             0.00000000e+00, 2.54716270e+04, -4.60117608e-01],
    #                             [3.38684249e+00, 3.47498246e-03, -6.35469633e-06, 6.96858127e-09,
    #                             -2.50658847e-12, -3.02081133e+04, 2.59023285e+00],
    #                             [3.21293640e+00, 1.12748635e-03, -5.75615047e-07, 1.31387723e-09,
    #                             -8.76855392e-13, -1.00524902e+03, 6.03473759e+00]])
    expected_nasa = {'O2': numpy.array([3.21293640e+00, 1.12748635e-03,
                                       -5.75615047e-07, 1.31387723e-09,
                                       -8.76855392e-13, -1.00524902e+03, 6.03473759e+00]),
                    'H2O': numpy.array([3.38684249e+00, 3.47498246e-03,
                                       -6.35469633e-06, 6.96858127e-09,
                                       -2.50658847e-12, -3.02081133e+04, 2.59023285e+00]),
                    'H': numpy.array([2.50000000e+00, 0.00000000e+00,
                                     0.00000000e+00, 0.00000000e+00,
                                     0.00000000e+00, 2.54716270e+04, -4.60117608e-01])}
    rev_rxn_obj = parser.reaction_list[0]
    #assert numpy.isclose(rev_rxn_obj.NASA_poly_coefs, expected_nasa).all()
    assert (numpy.isclose(rev_rxn_obj.NASA_poly_coefs['H2O'], expected_nasa['H2O'])).all()
    assert (numpy.isclose(rev_rxn_obj.NASA_poly_coefs['O2'], expected_nasa['O2'])).all()
    assert (numpy.isclose(rev_rxn_obj.NASA_poly_coefs['H'], expected_nasa['H'])).all()









# ======================= TESTS FOR REACTION OBJECT ====================== #

@pytest.fixture
def test_base_reaction():
    """Returns a valid reaction (from rxns.xml)"""
    return Reaction(rxn_type="Elementary",
                    is_reversible=False,
                    rxn_equation="H2 + OH =] H2O + H",
                    species_list=['H', 'O', 'OH', 'H2', 'H2O', 'O2'],
                    rate_coeffs_components={'k': 10},
                    reactant_stoich_coeffs={'H2' :1, 'OH':1},
                    product_stoich_coeffs={'H2O' :1, 'H':1})

def test_Reaction_special_fcns(test_base_reaction):
    """Test special functions for Reaction object"""
    info_reaction = "Reaction : H2 + OH =] H2O + H"
    assert str(test_base_reaction) == info_reaction
    assert len(test_base_reaction) == 4

def test_Reaction_get_unique_species(test_base_reaction):
    """Test get_unique_species function (used to get
    values of species_list)."""
    assert 'H2O' in test_base_reaction.unique_species
    assert 'H2' in test_base_reaction.unique_species
    assert 'H' in test_base_reaction.unique_species
    assert 'OH' in test_base_reaction.unique_species

def test_Reaction_set_zero_temperature(test_base_reaction):
    """Test setting reaction temperature to absolute 0."""
    with pytest.raises(ValueError):
        test_base_reaction.set_temperature(0)

def test_Reaction_set_neg_temperature(test_base_reaction):
    """Test setting reaction temperature to a negative value"""
    with pytest.raises(ValueError):
        test_base_reaction.set_temperature(-100)

def test_Reaction_set_valid_temperature(test_base_reaction):
    """Test setting reaction temperature to absolute 0."""
    test_base_reaction.set_temperature(10)
    assert test_base_reaction.temperature == 10

def test_Reaction_order_dictionaries(test_base_reaction):
    """Test ordering of reaction"""
    # order of species_list : ['H', 'O', 'OH', 'H2', 'H2O', 'O2']
    some_dict = {'H2':1, 'OH':2, 'H2O':3, 'H':4}
    expected = [4, 2, 1, 3]
    test_list = test_base_reaction.order_dictionaries(some_dict)
    assert test_list == expected

def test_Reaction_set_concentrations(test_base_reaction):
    """Test setting reaction with valid concentrations"""
    expected = [4, 2, 1, 3]
    test_base_reaction.set_concentrations({'H2':1, 'OH':2, 'H2O':3, 'H':4})
    assert (test_base_reaction.concentrations == expected).all()

def test_Reaction_set_neg_concentrations(test_base_reaction):
    """Test setting reaction with negative concentrations"""
    with pytest.raises(ValueError):
        test_base_reaction.set_concentrations({'H2':1, 'OH':2, 'H2O':3, 'H':-4})

def test_Reaction_compute_reaction_rate_coeff(test_base_reaction):
    """Test computing reaction rate coeff with bare Reaction."""
    with pytest.raises(NotImplementedError):
        test_base_reaction.compute_reaction_rate_coeff()

def test_Reaction_compute_progress_rate(test_base_reaction):
    """Test computing progress rate with bare Reaction."""
    with pytest.raises(NotImplementedError):
        test_base_reaction.compute_progress_rate()

def test_Reaction_compute_reaction_rate(test_base_reaction):
    """Test computing reaction rate with bare Reaction."""
    with pytest.raises(ReactionError):
        test_base_reaction.compute_reaction_rate()








# ======================= TESTS FOR IRREVERSIBLE REACTION ====================== #

@pytest.fixture
def test_irrev_reaction():
    """Returns a valid reaction (from rxns.xml)"""
    return IrreversibleReaction(rxn_type="Elementary",
                                is_reversible=False,
                                rxn_equation="H2 + OH =] H2O + H",
                                species_list=['H', 'O', 'OH', 'H2', 'H2O', 'O2'],
                                rate_coeffs_components={'k': 10},
                                reactant_stoich_coeffs={'H2' :1, 'OH':1},
                                product_stoich_coeffs={'H2O' :1, 'H':1})

@pytest.fixture
def test_irrev_reaction_arrhenius():
    """Returns a reaction with Arrhenius reaction rate coefficient."""
    return IrreversibleReaction(rxn_type="Elementary",
                                is_reversible=False,
                                rxn_equation="H2 + OH =] H2O + H",
                                species_list=['H', 'O', 'OH', 'H2', 'H2O', 'O2'],
                                rate_coeffs_components={'A': 10, 'E': 100},
                                reactant_stoich_coeffs={'H2' :1, 'OH':1},
                                product_stoich_coeffs={'H2O' :1, 'H':1})

@pytest.fixture
def test_irrev_reaction_modified_arr():
    """Returns a reaction with modified Arrhenius reaction rate coefficient."""
    return IrreversibleReaction(rxn_type="Elementary",
                                is_reversible=False,
                                rxn_equation="H2 + OH =] H2O + H",
                                species_list=['H', 'O', 'OH', 'H2', 'H2O', 'O2'],
                                rate_coeffs_components={'A': 10, 'E': 100, 'b':0.5},
                                reactant_stoich_coeffs={'H2' :1, 'OH':1},
                                product_stoich_coeffs={'H2O' :1, 'H':1})

def test_wrongly_classified_irrev_rxn_reversible():
    """Test initializing wrongly classified IrreversibleReaction (actually reversible)"""
    with pytest.raises(IrreversibleReactionError):
        test = IrreversibleReaction(rxn_type="Elementary",
                                    is_reversible=True,
                                    rxn_equation="H2 + OH =] H2O + H",
                                    species_list=['H', 'O', 'OH', 'H2', 'H2O', 'O2'],
                                    rate_coeffs_components={'k': 10},
                                    reactant_stoich_coeffs={'H2' :1, 'OH':1},
                                    product_stoich_coeffs={'H2O' :1, 'H':1})

def test_wrongly_classified_irrev_rxn_nonElementary():
    """Test initializing wrongly classified IrreversibleReaction (actually non-elementary)"""
    with pytest.raises(IrreversibleReactionError):
        test = IrreversibleReaction(rxn_type="Non-elementary",
                                    is_reversible=False,
                                    rxn_equation="H2 + OH =] H2O + H",
                                    species_list=['H', 'O', 'OH', 'H2', 'H2O', 'O2'],
                                    rate_coeffs_components={'k': 10},
                                    reactant_stoich_coeffs={'H2' :1, 'OH':1},
                                    product_stoich_coeffs={'H2O' :1, 'H':1})

def test_IrrevReaction_no_set_concentrations(test_irrev_reaction):
    """Test computing progress rate without setting species concentrations."""
    with pytest.raises(ValueError):
        test_irrev_reaction.compute_progress_rate()

def test_IrrevReaction_compute_reaction_rate_coeff_constant(test_irrev_reaction):
    """Test reaction constant reaction rate coefficient"""
    k = test_irrev_reaction.compute_reaction_rate_coeff()
    assert k == 10

def test_IrrevReaction_compute_reaction_rate_coeff_invalid_constant():
    """Test reaction constant reaction rate coefficient
    but invalid constant (non-positive)"""
    test_rxn = IrreversibleReaction(rxn_type="Elementary",
                                    is_reversible=False,
                                    rxn_equation="H2 + OH =] H2O + H",
                                    species_list=['H', 'O', 'OH', 'H2', 'H2O', 'O2'],
                                    rate_coeffs_components={'k': -10},
                                    reactant_stoich_coeffs={'H2' :1, 'OH':1},
                                    product_stoich_coeffs={'H2O' :1, 'H':1})
    with pytest.raises(ValueError):
        test_rxn.compute_reaction_rate_coeff()

def test_IrrevReaction_compute_reaction_rate_coeff_arrhenius(test_irrev_reaction_arrhenius):
    """Test irrev reaction with arrhenius rate coefficient."""
    # when T not inputed by user
    try:
        k = test_irrev_reaction_arrhenius.compute_reaction_rate_coeff() 
    except ValueError:
        test_irrev_reaction_arrhenius.set_temperature(10)
        k = test_irrev_reaction_arrhenius.compute_reaction_rate_coeff()
    assert numpy.isclose(k, 3.0035490889639616)

def test_IrrevReaction_compute_reaction_rate_coeff_mod_arrhenius(test_irrev_reaction_modified_arr):
    """Test reaction with modified arrhenius rate coefficient."""
    # when T not inputed by user
    try:
        k = test_irrev_reaction_modified_arr.compute_reaction_rate_coeff() 
    except ValueError:
        test_irrev_reaction_modified_arr.set_temperature(10)
        k = test_irrev_reaction_modified_arr.compute_reaction_rate_coeff()
    assert numpy.isclose(k, 9.4980561852498244)

def test_IrrevReaction_compute_progress_rate():
    """Test compute_progress_rate() for an elementary,
    irreversible reaction."""
    test = IrreversibleReaction(rxn_type="Elementary",
                                is_reversible=False,
                                rxn_equation="A + B =] C",
                                species_list=['A', 'B', 'C'],
                                rate_coeffs_components={'k': 10},
                                reactant_stoich_coeffs={'A' :2, 'B':1},
                                product_stoich_coeffs={'C': 1})
    test.set_concentrations({'A': 1, 'B':2, 'C':3})
    w = test.compute_progress_rate()
    expected = numpy.array([ 20.])
    assert w == expected


def test_IrrevReaction_compute_reaction_rate(test_irrev_reaction):
    """Test compute_reaction_rate() for an elementary,
    irreversible reaction."""
    test = IrreversibleReaction(rxn_type="Elementary",
                                is_reversible=False,
                                rxn_equation="A + B =] C",
                                species_list=['A', 'B', 'C'],
                                rate_coeffs_components={'k': 10},
                                reactant_stoich_coeffs={'A' :2, 'B':1},
                                product_stoich_coeffs={'C': 1})
    test.set_concentrations({'A': 1, 'B':2, 'C':3})
    rxnrate = test.compute_reaction_rate()
    expected = numpy.array([-40.0, -20.0, 20.0])
    assert (rxnrate == expected).all()

def test_IrrevReaction_compute_reaction_rate_neg_reactant_stoich_coeffs(test_irrev_reaction):
    """Test compute_reaction_rate() for an elementary, irreversible reaction."""
    test = IrreversibleReaction(rxn_type="Elementary",
                                is_reversible=False,
                                rxn_equation="A + B =] C",
                                species_list=['A', 'B', 'C'],
                                rate_coeffs_components={'k': 10},
                                reactant_stoich_coeffs={'A' :2, 'B':-1},
                                product_stoich_coeffs={'C': 1})
    test.set_concentrations({'A': 1, 'B':2, 'C':3})
    with pytest.raises(ValueError):
        rxnrate = test.compute_reaction_rate()
    
def test_IrrevReaction_compute_reaction_rate_neg_product_stoich_coeffs(test_base_reaction):
    """Test compute_reaction_rate() for an elementary, irreversible reaction."""
    test = IrreversibleReaction(rxn_type="Elementary",
                                is_reversible=False,
                                rxn_equation="A + B =] C",
                                species_list=['A', 'B', 'C'],
                                rate_coeffs_components={'k': 10},
                                reactant_stoich_coeffs={'A' :2, 'B':1},
                                product_stoich_coeffs={'C': -1})
    test.set_concentrations({'A': 1, 'B':2, 'C':3})
    with pytest.raises(ValueError):
        rxnrate = test.compute_reaction_rate()












# ======================= TESTS FOR REVERSIBLE REACTION ====================== #

@pytest.fixture
def test_rev_reaction():
    """Returns a valid reaction (from rev_rxn.xml)"""
    return ReversibleReaction(rxn_type="Elementary",
                    is_reversible=True,
                    rxn_equation="H + O2 [=] H2O ",
                    species_list=['H', 'H2O', 'O2'],
                    rate_coeffs_components={'k': 10},
                    reactant_stoich_coeffs={'H' :2, 'O2':1},
                    product_stoich_coeffs={'H2O' :2})

def test_wrongly_classified_rev_rxn_irreversible():
    """Tests initializing wrongly classified ReversibleReaction (actually irreversible)"""
    with pytest.raises(ReversibleReactionError):
        test = ReversibleReaction(rxn_type="Elementary",
                                    is_reversible=False,
                                    rxn_equation="H2 + OH =] H2O + H",
                                    species_list=['H', 'O', 'OH', 'H2', 'H2O', 'O2'],
                                    rate_coeffs_components={'k': 10},
                                    reactant_stoich_coeffs={'H2' :1, 'OH':1},
                                    product_stoich_coeffs={'H2O' :1, 'H':1})

def test_wrongly_classified_rev_rxn_nonElementary():
    """Tests initializing wrongly classified ReversibleReaction (actually non-elementary)"""
    with pytest.raises(ReversibleReactionError):
        test = ReversibleReaction(rxn_type="Non-elementary",
                                    is_reversible=True,
                                    rxn_equation="H2 + OH =] H2O + H",
                                    species_list=['H', 'O', 'OH', 'H2', 'H2O', 'O2'],
                                    rate_coeffs_components={'k': 10},
                                    reactant_stoich_coeffs={'H2' :1, 'OH':1},
                                    product_stoich_coeffs={'H2O' :1, 'H':1})

def test_compute_reaction_rate_coeff_rev_rxn_without_setting_nasa(test_rev_reaction):
    """Tests trying to compute reaction rate coeffs without setting nasa coeffs"""
    with pytest.raises(ValueError):
        test_rev_reaction.compute_reaction_rate_coeff()


def test_compute_reaction_rate_coeff_rev_rxn(test_rev_reaction):
    """Tests computing reaction rate coeff for a reversible reaction"""
    T = 400
    test_rev_reaction.set_temperature(T)
    lowT_nasa = numpy.array([[2.50000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
                            0.00000000e+00, 2.54716270e+04, -4.60117608e-01],
                            [3.38684249e+00, 3.47498246e-03, -6.35469633e-06, 6.96858127e-09,
                            -2.50658847e-12, -3.02081133e+04, 2.59023285e+00],
                            [3.21293640e+00, 1.12748635e-03, -5.75615047e-07, 1.31387723e-09,
                            -8.76855392e-13, -1.00524902e+03, 6.03473759e+00]])
    test_rev_reaction.set_NASA_poly_coefs(lowT_nasa)
    kf, kb = test_rev_reaction.compute_reaction_rate_coeff(T=T)
    assert kf == 10
    assert numpy.isclose(kb, 0)

def test_compute_progress_rate_rev_rxn(test_rev_reaction):
    """Tests computing reaction rate coeff for a reversible reaction"""
    T = 400
    test_rev_reaction.set_temperature(T)
    lowT_nasa = numpy.array([[2.50000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
                            0.00000000e+00, 2.54716270e+04, -4.60117608e-01],
                            [3.38684249e+00, 3.47498246e-03, -6.35469633e-06, 6.96858127e-09,
                            -2.50658847e-12, -3.02081133e+04, 2.59023285e+00],
                            [3.21293640e+00, 1.12748635e-03, -5.75615047e-07, 1.31387723e-09,
                            -8.76855392e-13, -1.00524902e+03, 6.03473759e+00]])
    test_rev_reaction.set_NASA_poly_coefs(lowT_nasa)
    test_rev_reaction.set_concentrations(X={'H':1, 'O2':1, 'H2O':1})
    prog_rate = test_rev_reaction.compute_progress_rate(T=T)
    assert prog_rate == 10

def test_compute_progress_rate_rev_rxn_without_setting_concen(test_rev_reaction):
    """Tests computing reaction rate coeff for a reversible reaction
    WITHOUT setting concentration"""
    T = 400
    test_rev_reaction.set_temperature(T)
    lowT_nasa = numpy.array([[2.50000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
                            0.00000000e+00, 2.54716270e+04, -4.60117608e-01],
                            [3.38684249e+00, 3.47498246e-03, -6.35469633e-06, 6.96858127e-09,
                            -2.50658847e-12, -3.02081133e+04, 2.59023285e+00],
                            [3.21293640e+00, 1.12748635e-03, -5.75615047e-07, 1.31387723e-09,
                            -8.76855392e-13, -1.00524902e+03, 6.03473759e+00]])
    test_rev_reaction.set_NASA_poly_coefs(lowT_nasa)
    with pytest.raises(ValueError):
        prog_rate = test_rev_reaction.compute_progress_rate(T=T)




# ======================= TESTS FOR REACTIONCOEFF ====================== #


def test_ReactionCoeff_constant():
    """Test when reaction rate coefficient is constant"""
    k_parameters = {'k': 10}
    k_test = ReactionCoeff(k_parameters).k
    assert k_test == 10


def test_ReactionCoeff_constant_with_others():
    """Test when reaction rate coefficient and A', 'b', 'E' all exists"""
    k_parameters = {'k': 10,'A': 10**7, "q":8}
    with pytest.raises(ValueError) as excinfo:
        k_test = ReactionCoeff(k_parameters).k
        #assert excinfo.value.message == "Invalid key in the input!"

def test_ReactionCoeff_invalid_constant():
    """Test when reaction rate coefficient is constant but invalid (non-positive)"""
    k_parameters = {'k': -10}
    with pytest.raises(ValueError):
        k_test = ReactionCoeff(k_parameters).k

def test_ReactionCoeff_constant_with_T():
    """Test when reaction rate coefficient is constant but T entered for some reason
    (should have no effect)"""
    k_parameters = {'k': 10}
    T = 10
    k_test = ReactionCoeff(k_parameters, T).k
    assert k_test == 10

def test_ReactionCoeff_arrhenius():
    """Test when reaction rate coefficient is Arrhenius"""
    k_parameters = {'A': 10**7, 'E':10**3}
    T = 10**2
    k_test = ReactionCoeff(k_parameters, T).k
    assert numpy.isclose(k_test, 3003549.08896)

def test_ReactionCoeff_arrhenius_invalid_A():
    """Test when reaction rate coefficient is Arrhenius but A is invalid (non-positive)"""
    k_parameters = {'A': 0, 'E':100}
    T = 10
    with pytest.raises(ValueError):
        k_test = ReactionCoeff(k_parameters, T).k
   
def test_ReactionCoeff_arrhenius_T_not_set():
    """Test when reaction rate coefficient is Arrhenius but T is not set by user"""
    k_parameters = {'A': 10, 'E':100}
    with pytest.raises(ValueError):
        k_test = ReactionCoeff(k_parameters).k

def test_ReactionCoeff_arrhenius_invalid_T():
    """Test when reaction rate coefficient is Arrhenius but T is invalid (non-positive)"""
    k_parameters = {'A': 10, 'E':100}
    T = -10
    with pytest.raises(ValueError):
        k_test = ReactionCoeff(k_parameters, T).k

def test_ReactionCoeff_arrhenius_invalid_R():
    """Test when reaction rate coefficient is Arrhenius but R is invalid (non-positive)"""
    k_parameters = {'A': 10, 'E':100, 'R':-100}
    T = 10
    with pytest.raises(ValueError):
        k_test = ReactionCoeff(k_parameters, T).k

def test_ReactionCoeff_arrhenius_changing_R():
    """Test when reaction rate coefficient is Arrhenius but R is changed by user"""
    k_parameters = {'A': 10, 'E':100, 'R':10.453}
    T = 10
    with pytest.raises(UserWarning):
        k_test = ReactionCoeff(k_parameters, T).k

def test_ReactionCoeff_mod_arrhenius():
    """Test when reaction rate coefficient is modified Arrhenius"""
    k_parameters = {'A': 10**7, 'E':10**3, 'b': 0.5}
    T = 10**2
    k_test = ReactionCoeff(k_parameters, T).k
    assert numpy.isclose(k_test, 30035490.8896)

def test_ReactionCoeff_mod_arrhenius_invalid_A():
    """Test when reaction rate coefficient is modified
    Arrhenius but A is invalid (non-positive)"""
    k_parameters = {'A': -10, 'E':100, 'b':0.5}
    T = 10
    with pytest.raises(ValueError):
        k_test = ReactionCoeff(k_parameters, T).k

def test_ReactionCoeff_mod_arrhenius_invalid_b():
    """Test when reaction rate coefficient is modified
    Arrhenius but B is invalid (not real)"""
    k_parameters = {'A': 10, 'E':100, 'b':0.5j}
    T = 10
    with pytest.raises(TypeError):
        k_test = ReactionCoeff(k_parameters, T).k
   
def test_ReactionCoeff_mod_arrhenius_T_not_set():
    """Test when reaction rate coefficient is modified
    Arrhenius but T is not set by user"""
    k_parameters = {'A': 10, 'E':100, 'b':0.5}
    with pytest.raises(ValueError):
        k_test = ReactionCoeff(k_parameters).k

def test_ReactionCoeff_mod_arrhenius_invalid_T():
    """Test when reaction rate coefficient is modified
    Arrhenius but T is invalid (non-positive)"""
    k_parameters = {'A': 10, 'E':100, 'b':0.5}
    T = -10
    with pytest.raises(ValueError):
        k_test = ReactionCoeff(k_parameters, T).k

def test_ReactionCoeff_mod_arrhenius_invalid_R():
    """Test when reaction rate coefficient is modified 
    Arrhenius but R is invalid (non-positive)"""
    k_parameters = {'A': 10, 'E':100, 'R':-100, 'b':0.5}
    T = 10
    with pytest.raises(ValueError):
        k_test = ReactionCoeff(k_parameters, T).k

def test_ReactionCoeff_mod_arrhenius_changing_R():
    """Test when reaction rate coefficient is modified
    Arrhenius but R is changed by user"""
    k_parameters = {'A': 10, 'E':100, 'R':10.453, 'b':0.5}
    T = 10
    with pytest.raises(UserWarning):
        k_test = ReactionCoeff(k_parameters, T).k







# ======================= TESTS FOR BACKWARDCOEFF ====================== #

@pytest.fixture
def test_backward_coeff():
    """Returns a working (but artificial) example of
    backward reaction rate coefficient."""
    expected_nasa = numpy.array([[1,0,0,0,0,0,0],
                                [1,0,0,0,0,0,0],
                                [1,0,0,0,0,0,0]])
    k_f = 100
    nu_i = numpy.array([-2, -1, 2])
    bkwd_coeff = BackwardCoeff(nu_i, expected_nasa)
    return bkwd_coeff

def test_backwardCoeff_gamma(test_backward_coeff):
    """Tests value of gamma for working example."""
    assert test_backward_coeff.gamma == -1

def test_backwardCoeff_computing_H(test_backward_coeff):
    """Tests computing H/RT for working example."""
    T = 100
    expected_H_over_RT = numpy.array([1, 1, 1])
    assert numpy.isclose(test_backward_coeff.H_over_RT(T),
                         expected_H_over_RT).all()

def test_backwardCoeff_computing_H_neg_T(test_backward_coeff):
    """Tests computing H/RT for working example with neg T."""
    T = -100
    with pytest.raises(ValueError):
        test_backward_coeff.H_over_RT(T)

def test_backwardCoeff_computing_S(test_backward_coeff):
    """Tests computing S/R for working example."""
    T = 100
    expected_S_over_R = numpy.array([4.60517, 4.60517, 4.60517])
    assert numpy.isclose(test_backward_coeff.S_over_R(T),
                         expected_S_over_R).all()

def test_backwardCoeff_computing_S_neg_T(test_backward_coeff):
    """Tests computing S/R for working example with neg T."""
    T = -100
    with pytest.raises(ValueError):
        test_backward_coeff.S_over_R(T)

def test_backwardCoeff_computeCoeff(test_backward_coeff):
    """Tests computing k_b for working example."""
    T = 100
    k_f = 100
    expected_delta_S_over_R = -4.60517
    expected_delta_H_over_RT = -1

    fact =  test_backward_coeff.p0 / test_backward_coeff.R / T
    expected_gamma = -1
    expected_ke = (fact ** expected_gamma) * (numpy.exp(expected_delta_S_over_R - 
                                                        expected_delta_H_over_RT))

    expected_kb_val = 442457 # 100 / 2.260104919e-6

    assert numpy.isclose(test_backward_coeff.compute_backward_coeffs(k_f, T),
                         expected_kb_val)

def test_backwardCoeff_computeCoeff_neg_T(test_backward_coeff):
    """Tests computing k_b for working example with neg T."""
    T = -100
    k_f = 100
    with pytest.raises(ValueError):
        test_backward_coeff.compute_backward_coeffs(k_f, T)

















# def test_overall_workflow_elementary_rxn():
#     """Test overall workflow (by checking final reaction rate)
#     for an irreversible elementary reaction."""
#     xml_filename = "rxns.xml"
#     parser = ReactionParser(xml_filename)
#     parser()
#     rxn1 = parser.reaction_list[0]
#     rxn1.set_concentrations({'H':1, 'O2':2, 'OH':0, 'O':0, 'H2O':0, 'H2':0})
#     rxn1.set_temperature(100)
#     rxnrate = rxn1.compute_reaction_rate()
#     expected = 0.0
#     assert rxnrate == expected


# # TODO: Add test using xml file with a set of reactions that David gave us!



