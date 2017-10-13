#Generic Function Written for both Reaction/Progress Rate for System of one or more equations. This file consists
#the unit tests for the algorithm.

import chemkin as original
import numpy as np

def test_single():
    #Tests whether the generic algorithm works for a system of one equation
    assert original.prog_rate([[2.0,1.0,0.0]],None,[1.0,2.0,3.0],[10],False)==[20.0]
    
def test_mulitresult():
    #Tests whether the generic algorithm works for a system of multiple equations
    assert original.prog_rate([[1.0,2.0,0.0],[2.0,0.0,2.0]],[[0.0,0.0,2.0],[0.0,1.0,1.0]],[1.0,2.0,1.0],[10,10],False)==[40.0, 10.0]
    
def test_negative_k_error():
    #Tests whether the generic algorithm works for a system of multiple equations
    try:
        assert original.prog_rate([[2.0,1.0,0.0]],None,[1.0,2.0,3.0],[-10],False)
    except ValueError as err:
        assert(type(err) == ValueError)

def test_negative_mole_error():
    #Tests whether the generic algorithm works for a system of multiple equations
    try:
        assert original.prog_rate([[-2.0,1.0,0.0]],None,[1.0,2.0,3.0],[10],False)
    except ValueError as err:
        assert(type(err) == ValueError)

def test_negative_concentration_error():
    #Tests whether the generic algorithm works for a system of multiple equations
    try:
        assert original.prog_rate([[2.0,1.0,0.0]],None,[1.0,-2.0,3.0],[10],False)
    except ValueError as err:
        assert(type(err) == ValueError)