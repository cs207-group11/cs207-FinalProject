
import numpy as np
def prog_rate(v1,v2,x,k,check):
    """
    Returns the progress (and reaction rate) of a system of elementary chemical reactions. 
    
    INPUTS
    =======
    v1 - No. of moles of the reactants of the system of reactions where each reaction is an independent vector.
    v2 - No. of moles of the products of the system of reactions where each reaction is an independent vector.
    x - Vector of concentrations of each unique element
    k - Reaction rate coefficient
    check - boolean value to check if progress rate to be returned or both progress rate and reaction rate. True --> Both.
    
    RETURNS
    ========
    prate - Progress rate of the system of reactions.
    reaction rates of unique reactants in the systems of reactions.
    
    EXAMPLES
    =========
    >>> prog_rate([[2.0,1.0,0.0]],None,[1.0,2.0,3.0],[10],False)
    [20.0]
    >>> prog_rate([[1.0,2.0,0.0],[2.0,0.0,2.0]],[[0.0,0.0,2.0],[0.0,1.0,1.0]],[1.0,2.0,1.0],[10,10],False)
    [40.0, 10.0]
    >>> prog_rate([[1.0,2.0,0.0],[0.0,0.0,2.0]],[[0.0,0.0,1.0],[1.0,2.0,0.0]],[1.0,2.0,1.0],[10,10],True)
    ([40.0, 10.0], array([-30., -60.,  20.]))
    """
    for i in k:
        if i<0:
            raise ValueError("Reaction Rate Coefficient must be positive.")
    for i in x:
        if i<0:
            raise ValueError("Concentrations of unique reactants/products must be positive.")
    for i in v1:
        for j in i:
            if j<0:
                raise ValueError("No. of moles of a reactant cannot be negative.")
    if v2!=None:
        for i in v2:
            for j in i:
                if j<0:
                    raise ValueError("No. of moles of a Product cannot be negative.")
    if isinstance(check,bool)!=True:
        raise ValueError("Check must be boolean.")
        
    prate = k
    for i in range(len(v1)):
        reaction = v1[i]
        for index,value in enumerate(x):
            prate[i] *= value**(reaction[index])
    if check == False: 
        return prate
    else:
        return prate, np.dot((np.array(v2).T-np.array(v1).T),np.array(prate).T)

if __name__ == "__main__":
    import doctest
    doctest.testmod()
    
#The Unit tests for the generic algorithm are available below.