
"""Test module for Reactions."""

import numpy
import pytest
import warnings

from chemkinlib.reactions import Reactions, ReactionRateCoeffs

# Treat warnings like errors (for testing purposes)
warnings.simplefilter("error")


# ======================= TESTS FOR REACTION OBJECT ====================== #

@pytest.fixture
def test_base_reaction():
    """Returns a valid reaction (from rxns.xml)"""
    return Reactions.Reaction(rxn_type="Elementary",
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

def test_Reaction_set_concentrations_invalid(test_base_reaction):
    """Test setting reaction with invalid concentrations"""
    with pytest.raises(KeyError):
        test_base_reaction.set_concentrations({'H2':1, 'CH4':2, 'H2O':3, 'H':4})

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
    with pytest.raises(Reactions.ReactionError):
        test_base_reaction.compute_reaction_rate()








# ======================= TESTS FOR IRREVERSIBLE REACTION ====================== #

@pytest.fixture
def test_irrev_reaction():
    """Returns a valid reaction (from rxns.xml)"""
    return Reactions.IrreversibleReaction(rxn_type="Elementary",
                                is_reversible=False,
                                rxn_equation="H2 + OH =] H2O + H",
                                species_list=['H', 'O', 'OH', 'H2', 'H2O', 'O2'],
                                rate_coeffs_components={'k': 10},
                                reactant_stoich_coeffs={'H2' :1, 'OH':1},
                                product_stoich_coeffs={'H2O' :1, 'H':1})

@pytest.fixture
def test_irrev_reaction_arrhenius():
    """Returns a reaction with Arrhenius reaction rate coefficient."""
    return Reactions.IrreversibleReaction(rxn_type="Elementary",
                                is_reversible=False,
                                rxn_equation="H2 + OH =] H2O + H",
                                species_list=['H', 'O', 'OH', 'H2', 'H2O', 'O2'],
                                rate_coeffs_components={'A': 10, 'E': 100},
                                reactant_stoich_coeffs={'H2' :1, 'OH':1},
                                product_stoich_coeffs={'H2O' :1, 'H':1})

@pytest.fixture
def test_irrev_reaction_modified_arr():
    """Returns a reaction with modified Arrhenius reaction rate coefficient."""
    return Reactions.IrreversibleReaction(rxn_type="Elementary",
                                is_reversible=False,
                                rxn_equation="H2 + OH =] H2O + H",
                                species_list=['H', 'O', 'OH', 'H2', 'H2O', 'O2'],
                                rate_coeffs_components={'A': 10, 'E': 100, 'b':0.5},
                                reactant_stoich_coeffs={'H2' :1, 'OH':1},
                                product_stoich_coeffs={'H2O' :1, 'H':1})

def test_wrongly_classified_irrev_rxn_reversible():
    """Test initializing wrongly classified IrreversibleReaction (actually reversible)"""
    with pytest.raises(Reactions.IrreversibleReactionError):
        test = Reactions.IrreversibleReaction(rxn_type="Elementary",
                                    is_reversible=True,
                                    rxn_equation="H2 + OH =] H2O + H",
                                    species_list=['H', 'O', 'OH', 'H2', 'H2O', 'O2'],
                                    rate_coeffs_components={'k': 10},
                                    reactant_stoich_coeffs={'H2' :1, 'OH':1},
                                    product_stoich_coeffs={'H2O' :1, 'H':1})

def test_wrongly_classified_irrev_rxn_nonElementary():
    """Test initializing wrongly classified IrreversibleReaction (actually non-elementary)"""
    with pytest.raises(Reactions.IrreversibleReactionError):
        test = Reactions.IrreversibleReaction(rxn_type="Non-elementary",
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
    test_rxn = Reactions.IrreversibleReaction(rxn_type="Elementary",
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
    test = Reactions.IrreversibleReaction(rxn_type="Elementary",
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
    test = Reactions.IrreversibleReaction(rxn_type="Elementary",
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
    test = Reactions.IrreversibleReaction(rxn_type="Elementary",
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
    test = Reactions.IrreversibleReaction(rxn_type="Elementary",
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
    return Reactions.ReversibleReaction(rxn_type="Elementary",
                    is_reversible=True,
                    rxn_equation="H + O2 [=] H2O ",
                    species_list=['H', 'H2O', 'O2'],
                    rate_coeffs_components={'k': 10},
                    reactant_stoich_coeffs={'H' :2, 'O2':1},
                    product_stoich_coeffs={'H2O' :2})

def test_wrongly_classified_rev_rxn_irreversible():
    """Tests initializing wrongly classified ReversibleReaction (actually irreversible)"""
    with pytest.raises(Reactions.ReversibleReactionError):
        test = Reactions.ReversibleReaction(rxn_type="Elementary",
                                    is_reversible=False,
                                    rxn_equation="H2 + OH =] H2O + H",
                                    species_list=['H', 'O', 'OH', 'H2', 'H2O', 'O2'],
                                    rate_coeffs_components={'k': 10},
                                    reactant_stoich_coeffs={'H2' :1, 'OH':1},
                                    product_stoich_coeffs={'H2O' :1, 'H':1})

def test_wrongly_classified_rev_rxn_nonElementary():
    """Tests initializing wrongly classified ReversibleReaction (actually non-elementary)"""
    with pytest.raises(Reactions.ReversibleReactionError):
        test = Reactions.ReversibleReaction(rxn_type="Non-elementary",
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
    lowT_nasa = {'O2': numpy.array([3.21293640e+00, 1.12748635e-03,
                                       -5.75615047e-07, 1.31387723e-09,
                                       -8.76855392e-13, -1.00524902e+03, 6.03473759e+00]),
                    'H2O': numpy.array([3.38684249e+00, 3.47498246e-03,
                                       -6.35469633e-06, 6.96858127e-09,
                                       -2.50658847e-12, -3.02081133e+04, 2.59023285e+00]),
                    'H': numpy.array([2.50000000e+00, 0.00000000e+00,
                                     0.00000000e+00, 0.00000000e+00,
                                     0.00000000e+00, 2.54716270e+04, -4.60117608e-01])}
    # lowT_nasa = numpy.array([[2.50000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
    #                         0.00000000e+00, 2.54716270e+04, -4.60117608e-01],
    #                         [3.38684249e+00, 3.47498246e-03, -6.35469633e-06, 6.96858127e-09,
    #                         -2.50658847e-12, -3.02081133e+04, 2.59023285e+00],
    #                         [3.21293640e+00, 1.12748635e-03, -5.75615047e-07, 1.31387723e-09,
    #                         -8.76855392e-13, -1.00524902e+03, 6.03473759e+00]])
    test_rev_reaction.set_NASA_poly_coefs(lowT_nasa)
    kf, kb = test_rev_reaction.compute_reaction_rate_coeff(T=T)
    assert kf == 10
    assert numpy.isclose(kb, 0)

def test_compute_progress_rate_rev_rxn(test_rev_reaction):
    """Tests computing reaction rate coeff for a reversible reaction"""
    T = 400
    test_rev_reaction.set_temperature(T)
    # lowT_nasa = numpy.array([[2.50000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
    #                         0.00000000e+00, 2.54716270e+04, -4.60117608e-01],
    #                         [3.38684249e+00, 3.47498246e-03, -6.35469633e-06, 6.96858127e-09,
    #                         -2.50658847e-12, -3.02081133e+04, 2.59023285e+00],
    #                         [3.21293640e+00, 1.12748635e-03, -5.75615047e-07, 1.31387723e-09,
    #                         -8.76855392e-13, -1.00524902e+03, 6.03473759e+00]])
    lowT_nasa = {'O2': numpy.array([3.21293640e+00, 1.12748635e-03,
                                       -5.75615047e-07, 1.31387723e-09,
                                       -8.76855392e-13, -1.00524902e+03, 6.03473759e+00]),
                    'H2O': numpy.array([3.38684249e+00, 3.47498246e-03,
                                       -6.35469633e-06, 6.96858127e-09,
                                       -2.50658847e-12, -3.02081133e+04, 2.59023285e+00]),
                    'H': numpy.array([2.50000000e+00, 0.00000000e+00,
                                     0.00000000e+00, 0.00000000e+00,
                                     0.00000000e+00, 2.54716270e+04, -4.60117608e-01])}
    test_rev_reaction.set_NASA_poly_coefs(lowT_nasa)
    test_rev_reaction.set_concentrations(X={'H':1, 'O2':1, 'H2O':1})
    prog_rate = test_rev_reaction.compute_progress_rate(T=T)
    assert prog_rate == 10

def test_compute_progress_rate_rev_rxn_without_setting_concen(test_rev_reaction):
    """Tests computing reaction rate coeff for a reversible reaction
    WITHOUT setting concentration"""
    T = 400
    test_rev_reaction.set_temperature(T)
    # lowT_nasa = numpy.array([[2.50000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
    #                         0.00000000e+00, 2.54716270e+04, -4.60117608e-01],
    #                         [3.38684249e+00, 3.47498246e-03, -6.35469633e-06, 6.96858127e-09,
    #                         -2.50658847e-12, -3.02081133e+04, 2.59023285e+00],
    #                         [3.21293640e+00, 1.12748635e-03, -5.75615047e-07, 1.31387723e-09,
    #                         -8.76855392e-13, -1.00524902e+03, 6.03473759e+00]])
    lowT_nasa = {'O2': numpy.array([3.21293640e+00, 1.12748635e-03,
                                       -5.75615047e-07, 1.31387723e-09,
                                       -8.76855392e-13, -1.00524902e+03, 6.03473759e+00]),
                    'H2O': numpy.array([3.38684249e+00, 3.47498246e-03,
                                       -6.35469633e-06, 6.96858127e-09,
                                       -2.50658847e-12, -3.02081133e+04, 2.59023285e+00]),
                    'H': numpy.array([2.50000000e+00, 0.00000000e+00,
                                     0.00000000e+00, 0.00000000e+00,
                                     0.00000000e+00, 2.54716270e+04, -4.60117608e-01])}
    test_rev_reaction.set_NASA_poly_coefs(lowT_nasa)
    with pytest.raises(ValueError):
        prog_rate = test_rev_reaction.compute_progress_rate(T=T)
