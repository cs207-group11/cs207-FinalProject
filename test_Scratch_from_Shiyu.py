from Scratch_from_Shiyu import *

#test for correct results of computation
def test_const_result():
    coeff = ReactionCoeff()
    assert coeff.get_coeff(k=10, T=10e2) == 10


def test_Arrhenius_result():
    coeff = ReactionCoeff()
    assert (coeff.get_coeff(A=10e7, E=10e3, T=10e2) == 30035490.889639616)


def test_mod_Arrhenius_result():
    coeff = ReactionCoeff()
    assert (coeff.get_coeff(A=10e7, b=0.5, E=10e3, T=10e2) == 949805618.52498233)



#test for correct input
#arr(self,A, E, T, R=8.314)
def test_arr_value_A():
    coeff = ReactionCoeff()
    try:
        coeff.get_coeff(A=-1, E=10 ** 3, T=10 ** 2)
    except ValueError as err:
        assert (type(err) == ValueError)


def test_arr_value_T():
    coeff = ReactionCoeff()
    try:
        coeff.get_coeff(A=10 ** 7, E=10 ** 3, T=-1)
    except ValueError as err:
        assert (type(err) == ValueError)


def test_arr_value_R():
    coeff = ReactionCoeff()
    try:
        coeff.get_coeff(A=1, E=0.5, T=10 ** 3, R=-1)
    except ValueError as err:
        assert (type(err) == ValueError)

#mod_arr(self, A, b, E, T, R=8.314)
def test_mod_arr_value_A():
    coeff = ReactionCoeff()
    try:
        coeff.get_coeff(A=-1, b=0.5, E=10 ** 3, T=10 ** 2)
    except ValueError as err:
        assert (type(err) == ValueError)


def test_mod_arr_value_T():
    coeff = ReactionCoeff()
    try:
        coeff.get_coeff(A=10 ** 7, b=0.5, E=10 ** 3, T=-1)
    except ValueError as err:
        assert (type(err) == ValueError)


def test_mod_arr_value_R():
    coeff = ReactionCoeff()
    try:
        coeff.get_coeff(A=1, b=0.5, E=10 ** 3, T=10 ** 2, R=-1)
    except ValueError as err:
        assert (type(err) == ValueError)


def test_mod_arr_value_b():
    coeff = ReactionCoeff()
    try:
        coeff.get_coeff(10 ** 7, "1", 10 ** 3, 1)
    except TypeError as err:
        assert (type(err) == TypeError)