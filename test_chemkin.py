"""Test module for chemkin."""

import numpy
import pytest
import warnings
from chemkin import *

# Treat warnings like errors (for testing purposes)
warnings.simplefilter("error")


# # ======================= TESTS FOR REACTIONSYSTEM OBJECT ====================== #

# @pytest.fixture
# def test_base_reaction_system():
#     """Returns a valid reaction system"""
#     rxn1 = Reaction(rxn_type="Elementary",
#                     is_reversible=False,
#                     rxn_equation="H2 + OH =] H2O + H",
#                     species_list=['H', 'O', 'OH', 'H2', 'H2O', 'O2'],
#                     rate_coeffs_components={'k': 10},
#                     reactant_stoich_coeffs={'H2' :1, 'OH':1},
#                     product_stoich_coeffs={'H2O' :1, 'H':1})
#     rxn2 = Reaction(rxn_type="Elementary",
#                     is_reversible=True,
#                     rxn_equation="O2 + H =] H2O",
#                     species_list=['H', 'O', 'OH', 'H2', 'H2O', 'O2'],
#                     rate_coeffs_components={'k': 10},
#                     reactant_stoich_coeffs={'O2' :1, 'H':1},
#                     product_stoich_coeffs={'H2O' :1})
#     reaction_list = [rxn1, rxn2]
#     return ReactionSystem(reaction_list)

# def test_set_valid_temperature(test_base_reaction_system):
#     """Tests setting valid temperature to reaction system."""
#     test_sys = test_base_reaction_system
#     test_sys.set_temperature(100)
#     for rxnObj in test_sys.reaction_list:
#         assert rxnObj.temperature == 100

# def test_set_invalid_temperature(test_base_reaction_system):
#     """Tests setting valid temperature to reaction system."""
#     test_sys = test_base_reaction_system
#     with pytest.raises(ValueError):
#         test_sys.set_temperature(-100)

# def test_set_valid_concentrations(test_base_reaction_system):
#     """Tests setting valid concentrations to reaction system."""
#     test_sys = test_base_reaction_system
#     X = {'H': 1, 'OH': 1, 'H2': 1, 'O':1, 'H2O':1, 'O2': 1}
#     test_sys.set_concentrations(X)
#     for rxnObj in test_sys.reaction_list:
#         list_concentrations = rxnObj.order_dictionaries(X)
#         assert (list_concentrations == rxnObj.concentrations).all()

# def test_set_invalid_concentrations(test_base_reaction_system):
#     """Tests for when non-positive values for concentrations entered.""" 
#     test_sys = test_base_reaction_system
#     X = {'H': 1, 'OH': -1, 'H2': 1, 'O':1, 'H2O':1, 'O2': 1}
#     with pytest.raises(ValueError):
#         test_sys.set_concentrations(X)

# def test_set_incomplete_concentrations(test_base_reaction_system):
#     test_sys = test_base_reaction_system
#     X = {'H': 1, 'H2O': 1, 'O':1, 'H2O':1, 'O2': 1}
#     with pytest.raises(ValueError):
#         test_sys.set_concentrations(X)

# def test_set_wrong_concentrations(test_base_reaction_system):
#     """Tests for when wrong species (species not in reaction system)
#     assigned concentrations."""
#     test_sys = test_base_reaction_system
#     X = {'H': 1, 'H2O': 1, 'O':1, 'H2O':1, 'O2': 1, 'A': 1}
#     with pytest.raises(ValueError):
#         test_sys.set_concentrations(X)




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


def test_wrongly_classified_rev_rxn_irreversible():
    """Test initializing wrongly classified ReversibleReaction (actually irreversible)"""
    with pytest.raises(ReversibleReactionError):
        test = ReversibleReaction(rxn_type="Elementary",
                                    is_reversible=False,
                                    rxn_equation="H2 + OH =] H2O + H",
                                    species_list=['H', 'O', 'OH', 'H2', 'H2O', 'O2'],
                                    rate_coeffs_components={'k': 10},
                                    reactant_stoich_coeffs={'H2' :1, 'OH':1},
                                    product_stoich_coeffs={'H2O' :1, 'H':1})

def test_wrongly_classified_rev_rxn_nonElementary():
    """Test initializing wrongly classified ReversibleReaction (actually non-elementary)"""
    with pytest.raises(ReversibleReactionError):
        test = ReversibleReaction(rxn_type="Non-elementary",
                                    is_reversible=True,
                                    rxn_equation="H2 + OH =] H2O + H",
                                    species_list=['H', 'O', 'OH', 'H2', 'H2O', 'O2'],
                                    rate_coeffs_components={'k': 10},
                                    reactant_stoich_coeffs={'H2' :1, 'OH':1},
                                    product_stoich_coeffs={'H2O' :1, 'H':1})











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



