from chemkin import *
import numpy as np

def test_reaction_rate():
    v1 = np.array([[1,2,0],[0,0,2]])
    v2 = np.array([[0,0,1],[1,2,0]])
    x = np.array([1,2,1])
    k = np.array([10,10])
    assert(np.array_equal(reaction_rate(v1,v2,x,k), np.array([-30, -60, 20])))
    
    
def test_progress_rate():
    v = np.array([[1,2,0], [2, 0, 2]])
    x = np.array([1,2,1])
    k = [10,10]
    assert(progress_rate(v, x, k)[0] == 40 and progress_rate(v, x, k)[1] == 10)