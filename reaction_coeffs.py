
"""This module contains functions that return reaction rate coefficients
   
   MEMBERS:
   ========
   const:   Returns a constant reaction rate coefficient
   arr:     Returns the Arrhenius reaction rate coefficient
   mod_arr: Returns the modified Arrhenius reaction rate coefficient
"""
import numpy as np

def const(k=1.0):
    """Returns constant reaction rate coefficients k.
    
    INPUTS
    =======
    k: numerical value, no default value
    
    RETURNS
    ========
    k: constant reaction rate coefficients.
    
    EXAMPLES
    =========
    >>> const(10)
    10
    """
    return(k)

def arr(A,E,T,R=8.314):
    """Returns Arrhenius reaction rate coefficients k.
    
    INPUTS
    =======
    A: float, strictly positive, no default value
       The Arrhenius prefactor 
    E: float, no default value
       The Arrhenius parameter 
    T: float, strictly positive, no default value
       Temperature T, asuuming a Kelvin scale
    R: float, default value is 8.314, cannot be changed except to convert units
       The ideal gas constant 
    
    RETURNS
    ========
    k: Arrhenius reaction rate coefficients k,
       floats
       unless A or T is not postive
       in which case a ValueError exception is raised
    
    EXAMPLES
    =========
    >>> arr(10**7,10**3,10**2)
    3003549.0889639617
    """
    if (A <=0):
        raise ValueError("Arrhenius prefactor A must be postive!")

    if (T <=0):
        raise ValueError("Temperatures T must be postive!")
    
    if (R <=0):
        raise ValueError("Gas constant R must be postive!")
    
    k = A * np.exp(-E / R / T)
    return(k)
        
def mod_arr(A,b,E,T,R=8.314):
    """Returns Arrhenius reaction rate coefficients k.
    
    INPUTS
    =======
    A: float, strictly positive, no default value
       The Arrhenius prefactor 
    b: real, no default value
       Modified Arrhenius parameter
    E: float, no default value
       The Arrhenius parameter 
    T: float, strictly positive, no default value
       Temperature T, asuuming a Kelvin scale
    R: float, default value is 8.314, cannot be changed except to convert units
       The ideal gas constant 
    
    RETURNS
    ========
    k: Arrhenius reaction rate coefficients k,
       floats
       unless A or T is not postive
       in which case a ValueError exception is raised
       Or b is not a real number
       in which case a TypeError exception is raised
    
    EXAMPLES
    =========
    >>> mod_arr(10**7,0.5,10**3,10**2)
    30035490.889639616
    """
    import numbers
    
    if (A <= 0):
        raise ValueError("Parameter A must be postive!")
    if (T <= 0):
        raise ValueError("Parameter T must be postive!")         
    if (isinstance(b, numbers.Real)) == False:
        raise TypeError("Parameter b must be a real number!")  
    if (R <=0):
        raise ValueError("Gas constant R must be postive!")
        
    k = A * T**b * np.exp(-E / R / T)
    return(k)