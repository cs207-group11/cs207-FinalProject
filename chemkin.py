import numpy as np

def progress_rate(v, x, k):
    """Returns the progress rate for a system of reactions
   
    INPUTS
    =======
    v: The matrix of Stoichiometric coefficients of reactants
    x: The vector of concentration of species
    k: The vector of reaction rate coefficient
    
    RETURNS
    ========
    w: Progress rate of the reaction system
    
    EXAMPLES
    =========
    >>> progress_rate(np.array([[1,2,0], [2, 0, 2]]), np.array([1,2,1]), np.array([10,10]))
    np.array([ 40.,  10.])
    """
    w = np.zeros(len(v))
    for i in range(len(v)):
        w[i] = k[i] * np.prod(np.power(x, v[i]))
    return w

def reaction_rate(v1, v2, x, k):
    """Returns the progress rate for a system of reactions 
   
    INPUTS
    =======
    v1: The matrix of Stoichiometric coefficients of reactants
    v2: The matrix of Stoichiometric coefficients of products
    x: The vector of concentration of species
    k: The vector of reaction rate coefficient
    
    RETURNS
    ========
    f: Reaction rate of the reaction system
    
    EXAMPLES
    =========
    >>> reaction_rate(np.array([[1,2,0],[0,0,2]]), np.array([[0,0,1],[1,2,0]]), np.array([1,2,1]), np.array([10,10]))
    np.array([-30, -60, 20])
    """
    w = progress_rate(v1, x, k)
    f = np.zeros(len(x))
    for i in range(len(x)):
        f[i] = np.sum(v2[:,i]*w - v1[:,i]*w)
    return f