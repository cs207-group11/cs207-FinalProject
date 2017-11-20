
import numpy
import os
import pytest
import warnings

from chemkinlib.config import DATA_DIRECTORY
from chemkinlib.utils import Parser 
from chemkinlib.reactions import Reactions, ReactionSystems


# Treat warnings like errors (for testing purposes)
warnings.simplefilter("error")


@pytest.fixture
def test_rxn_sys():
    """Returns a valid reaction system"""
    xml_filename =  os.path.join(DATA_DIRECTORY, "rxn.xml")
    parser = Parser.ReactionParser(xml_filename)
    temp = 500 # "low" temperature range in NASA coeffs database
    concentrations = {'H':1, 'O2':1, 'H2O':1}
    rxnsys = ReactionSystems.ReactionSystem(parser.reaction_list, parser.NASA_poly_coefs, temp, concentrations)
    return rxnsys

def test_rxn_sys_invalid_temperature():
    """Tests setting up reaction system with invalid temperatures."""
    xml_filename =  os.path.join(DATA_DIRECTORY, "rxns_mixed.xml")
    parser = Parser.ReactionParser(xml_filename)
    concentrations = {'H':1, 'O2':2, 'OH':1, 'O':4, 'H2O':0, 'H2':1}
    temp = 0
    with pytest.raises(ValueError):
        rxnsys = ReactionSystems.ReactionSystem(parser.reaction_list, parser.NASA_poly_coefs, temp, concentrations)
    temp = -100
    with pytest.raises(ValueError):
        rxnsys = ReactionSystems.ReactionSystem(parser.reaction_list, parser.NASA_poly_coefs, temp, concentrations)


def test_rxn_sys_get_reaction_rate_for_1_rxn(test_rxn_sys):
    """Tests function to get reaction rate for a given system of reactions (just 1 reaction)."""
    # print(test_rxn_sys.involved_species)
    # assert test_rxn_sys.get_reaction_rate() == -15.
    rates = test_rxn_sys.sort_reaction_rates()
    assert rates['H'] == -30.
    assert rates['O2'] == -15.
    assert rates['H2O'] == 30.


def test_rxn_sys_get_reaction_rate_for_3_rxns():
    """Tests function to get reaction rate for a given system of reactions (more than 1 reaction)."""
    xml_filename =  os.path.join(DATA_DIRECTORY, "rxnsys.xml")
    parser = Parser.ReactionParser(xml_filename)
    temp = 10
    concentrations = {'H':1, 'O2':1, 'OH':1, 'O':1, 'H2O':1, 'H2':1}
    rxnsys = ReactionSystems.ReactionSystem(parser.reaction_list, parser.NASA_poly_coefs, temp, concentrations)
    rates = rxnsys.sort_reaction_rates()
    # assert rxnsys.get_reaction_rate() == -10.
    assert rates['H'] == -10.
    assert rates['O2'] == -15.
    assert rates['H2O'] == 40.
    assert rates['H2'] == -20.
    assert rates['O'] == -10.
    assert rates['OH'] == 0.


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
    xml_filename =  os.path.join(DATA_DIRECTORY, "rxn.xml")
    parser = Parser.ReactionParser(xml_filename)
    temp = 5000 # "high" temperature range in NASA coeffs database
    concentrations = {'H':1, 'O2':1, 'H2O':1}
    rxnsys = ReactionSystems.ReactionSystem(parser.reaction_list, parser.NASA_poly_coefs, temp, concentrations)
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
    xml_filename =  os.path.join(DATA_DIRECTORY, "rev_rxn.xml")
    parser = Parser.ReactionParser(xml_filename)
    temp = 500 # "low" temperature range in NASA coeffs database
    concentrations = {'H':1, 'O2':1, 'H2O':1}
    rxnsys = ReactionSystems.ReactionSystem(parser.reaction_list, parser.NASA_poly_coefs, temp, concentrations)
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
    assert (numpy.isclose(rev_rxn_obj.NASA_poly_coefs_dict['H2O'], expected_nasa['H2O'])).all()
    assert (numpy.isclose(rev_rxn_obj.NASA_poly_coefs_dict['O2'], expected_nasa['O2'])).all()
    assert (numpy.isclose(rev_rxn_obj.NASA_poly_coefs_dict['H'], expected_nasa['H'])).all()
