
"""Test module for ReactionRateCoeffs."""

import numpy
import pytest
import warnings

from chemkinlib.reactions import ReactionRateCoeffs 

# Treat warnings like errors (for testing purposes)
warnings.simplefilter("error")


# ======================= TESTS FOR REACTIONCOEFF ====================== #


def test_ReactionCoeff_constant():
    """Test when reaction rate coefficient is constant"""
    k_parameters = {'k': 10}
    k_test = ReactionRateCoeffs.ReactionCoeff(k_parameters).k
    assert k_test == 10


def test_ReactionCoeff_constant_with_others():
    """Test when reaction rate coefficient and A', 'b', 'E' all exists"""
    k_parameters = {'k': 10,'A': 10**7, "q":8}
    with pytest.raises(ValueError) as excinfo:
        k_test = ReactionRateCoeffs.ReactionCoeff(k_parameters).k
        #assert excinfo.value.message == "Invalid key in the input!"

def test_ReactionCoeff_invalid_constant():
    """Test when reaction rate coefficient is constant but invalid (non-positive)"""
    k_parameters = {'k': -10}
    with pytest.raises(ValueError):
        k_test = ReactionRateCoeffs.ReactionCoeff(k_parameters).k

def test_ReactionCoeff_constant_with_T():
    """Test when reaction rate coefficient is constant but T entered for some reason
    (should have no effect)"""
    k_parameters = {'k': 10}
    T = 10
    k_test = ReactionRateCoeffs.ReactionCoeff(k_parameters, T).k
    assert k_test == 10

def test_ReactionCoeff_arrhenius():
    """Test when reaction rate coefficient is Arrhenius"""
    k_parameters = {'A': 10**7, 'E':10**3}
    T = 10**2
    k_test = ReactionRateCoeffs.ReactionCoeff(k_parameters, T).k
    assert numpy.isclose(k_test, 3003549.08896)

def test_ReactionCoeff_arrhenius_invalid_A():
    """Test when reaction rate coefficient is Arrhenius but A is invalid (non-positive)"""
    k_parameters = {'A': 0, 'E':100}
    T = 10
    with pytest.raises(ValueError):
        k_test = ReactionRateCoeffs.ReactionCoeff(k_parameters, T).k
   
def test_ReactionCoeff_arrhenius_T_not_set():
    """Test when reaction rate coefficient is Arrhenius but T is not set by user"""
    k_parameters = {'A': 10, 'E':100}
    with pytest.raises(ValueError):
        k_test = ReactionRateCoeffs.ReactionCoeff(k_parameters).k

def test_ReactionCoeff_arrhenius_invalid_T():
    """Test when reaction rate coefficient is Arrhenius but T is invalid (non-positive)"""
    k_parameters = {'A': 10, 'E':100}
    T = -10
    with pytest.raises(ValueError):
        k_test = ReactionRateCoeffs.ReactionCoeff(k_parameters, T).k

def test_ReactionCoeff_arrhenius_invalid_R():
    """Test when reaction rate coefficient is Arrhenius but R is invalid (non-positive)"""
    k_parameters = {'A': 10, 'E':100, 'R':-100}
    T = 10
    with pytest.raises(ValueError):
        k_test = ReactionRateCoeffs.ReactionCoeff(k_parameters, T).k

def test_ReactionCoeff_arrhenius_changing_R():
    """Test when reaction rate coefficient is Arrhenius but R is changed by user"""
    k_parameters = {'A': 10, 'E':100, 'R':10.453}
    T = 10
    with pytest.raises(UserWarning):
        k_test = ReactionRateCoeffs.ReactionCoeff(k_parameters, T).k

def test_ReactionCoeff_mod_arrhenius():
    """Test when reaction rate coefficient is modified Arrhenius"""
    k_parameters = {'A': 10**7, 'E':10**3, 'b': 0.5}
    T = 10**2
    k_test = ReactionRateCoeffs.ReactionCoeff(k_parameters, T).k
    assert numpy.isclose(k_test, 30035490.8896)

def test_ReactionCoeff_mod_arrhenius_invalid_A():
    """Test when reaction rate coefficient is modified
    Arrhenius but A is invalid (non-positive)"""
    k_parameters = {'A': -10, 'E':100, 'b':0.5}
    T = 10
    with pytest.raises(ValueError):
        k_test = ReactionRateCoeffs.ReactionCoeff(k_parameters, T).k

def test_ReactionCoeff_mod_arrhenius_invalid_b():
    """Test when reaction rate coefficient is modified
    Arrhenius but B is invalid (not real)"""
    k_parameters = {'A': 10, 'E':100, 'b':0.5j}
    T = 10
    with pytest.raises(TypeError):
        k_test = ReactionRateCoeffs.ReactionCoeff(k_parameters, T).k
   
def test_ReactionCoeff_mod_arrhenius_T_not_set():
    """Test when reaction rate coefficient is modified
    Arrhenius but T is not set by user"""
    k_parameters = {'A': 10, 'E':100, 'b':0.5}
    with pytest.raises(ValueError):
        k_test = ReactionRateCoeffs.ReactionCoeff(k_parameters).k

def test_ReactionCoeff_mod_arrhenius_invalid_T():
    """Test when reaction rate coefficient is modified
    Arrhenius but T is invalid (non-positive)"""
    k_parameters = {'A': 10, 'E':100, 'b':0.5}
    T = -10
    with pytest.raises(ValueError):
        k_test = ReactionRateCoeffs.ReactionCoeff(k_parameters, T).k

def test_ReactionCoeff_mod_arrhenius_invalid_R():
    """Test when reaction rate coefficient is modified 
    Arrhenius but R is invalid (non-positive)"""
    k_parameters = {'A': 10, 'E':100, 'R':-100, 'b':0.5}
    T = 10
    with pytest.raises(ValueError):
        k_test = ReactionRateCoeffs.ReactionCoeff(k_parameters, T).k

def test_ReactionCoeff_mod_arrhenius_changing_R():
    """Test when reaction rate coefficient is modified
    Arrhenius but R is changed by user"""
    k_parameters = {'A': 10, 'E':100, 'R':10.453, 'b':0.5}
    T = 10
    with pytest.raises(UserWarning):
        k_test = ReactionRateCoeffs.ReactionCoeff(k_parameters, T).k







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
    bkwd_coeff = ReactionRateCoeffs.BackwardCoeff(nu_i, expected_nasa)
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
