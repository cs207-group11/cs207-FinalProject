import numpy as np
class ReactionRate():
    def __init__(self):
        raise NotImplementedError

    def progress_rate(self):
        raise NotImplementedError

    def reaction_rate(self):
        raise NotImplementedError


class IrreversibleReactionRate(ReactionRate):
    def __init__(self):
        pass

    def progress_rate(self, v, x, k):
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
        >>> r = IrreversibleReactionRate()
        >>> r.progress_rate(np.array([[1,2,0], [2, 0, 2]]), np.array([1,2,1]), np.array([10,10]))
        array([ 40.,  10.])
        """
        w = np.zeros(len(v))
        for i in range(len(v)):
            w[i] = k[i] * np.prod(np.power(x, v[i]))
        return w

    def reaction_rate(self, v1, v2, x, k):
        """Returns the progress rate for a system of reactions

        INPUTS
        =======
        v1: The matrix of Stoichiometric coefficients of reactants
        v2: The matrix of Stoichiometric coefficients of products
        x: The vector of concentration of species
        k: The vector of reaction rate coefficient

        RETURNS
        ========
        f: The vector of the reaction rate of each species in the reaction system

        EXAMPLES
        =========
        >>> r = IrreversibleReactionRate()
        >>> r.reaction_rate(np.array([[1,2,0],[0,0,2]]), np.array([[0,0,1],[1,2,0]]), np.array([1,2,1]), np.array([10,10]))
        array([-30., -60.,  20.])
        """
        w = self.progress_rate(v1, x, k)
        f = np.zeros(len(x))
        for i in range(len(x)):
            f[i] = np.sum(v2[:,i]*w - v1[:,i]*w)
        return f


class ReversibleReactionRate(ReactionRate):
    def progress_rate(self):
        raise NotImplementedError

    def reaction_rate(self):
        raise NotImplementedError


class ReactionCoeff():
    """This module contains functions that return reaction rate coefficients

       MEMBERS:
       ========
       const:   Returns a constant reaction rate coefficient
       arr:     Returns the Arrhenius reaction rate coefficient
       mod_arr: Returns the modified Arrhenius reaction rate coefficient
    """
    def __init__(self):
        pass

    def get_coeff(self, T, **kwargs):
        """
        :param kwargs: k, A, b, E ,T, R(optional, usually not changed)
        :return: reaction coefficient
        """
        if "k" in kwargs:
            return self.const(kwargs['k'])
        if "A" in kwargs and "E" in kwargs and "b" not in kwargs:
            if "R" in kwargs:
                return self.arr(A = kwargs['A'], E = kwargs['E'], T = T, R = kwargs['R'])
            else:
                return self.arr(A=kwargs['A'], E=kwargs['E'], T=T)
        if "A" in kwargs and "E" in kwargs  and "b" in kwargs:
            if "R" in kwargs:
                return self.mod_arr(A = kwargs['A'], E = kwargs['E'], R = kwargs['R'], b = kwargs['b'], T=T)
            else:
                return self.mod_arr(A=kwargs['A'], E=kwargs['E'], b=kwargs['b'], T=T)

    def const(self,k=1.0):
        """Returns constant reaction rate coefficients k.

        INPUTS
        =======
        k: numerical value, no default value

        RETURNS
        ========
        k: constant reaction rate coefficients.

        EXAMPLES
        =========
        >>> rc = ReactionCoeff()
        >>> rc.const(10)
        10
        """
        return (k)

    def arr(self,A, E, T, R=8.314):
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
        >>> rc = ReactionCoeff()
        >>> rc.arr(10**7,10**3,10**2)
        3003549.0889639617
        """
        if (A <= 0):
            raise ValueError("Arrhenius prefactor A must be postive!")

        if (T <= 0):
            raise ValueError("Temperatures T must be postive!")

        if (R <= 0):
            raise ValueError("Gas constant R must be postive!")

        k = A * np.exp(-E / R / T)
        return (k)

    def mod_arr(self, A, b, E, T, R=8.314):
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
        >>> rc = ReactionCoeff()
        >>> rc.mod_arr(10**7,0.5,10**3,10**2)
        30035490.889639616
        """
        import numbers

        if (A <= 0):
            raise ValueError("Parameter A must be postive!")
        if (T <= 0):
            raise ValueError("Parameter T must be postive!")
        if (isinstance(b, numbers.Real)) == False:
            raise TypeError("Parameter b must be a real number!")
        if (R <= 0):
            raise ValueError("Gas constant R must be postive!")

        k = A * T ** b * np.exp(-E / R / T)
        return (k)



##Simple test to see the whole structure works

rate_calculator = IrreversibleReactionRate()
coeff = ReactionCoeff()
T = [750, 1500, 2500]
v1 = np.array([[2,1,0,0,0],[0,0,1,1,0],[0,1,0,0,1]])
v2 = np.array([[1,0,2,0,0],[0,1,0,0,1],[0,0,1,1,0]])
x = np.array([2,1,0.5,1,1])
A1, b1, E1 = 1e8, 0.5, 5*1e4
k2 = coeff.const(1e4)
A3, E3 = 1e7, 1e4
reaction_rate = []
reaction_component = ["H2", "O2", "OH", "HO2", "H2O"]
for t in T:
    k1 = coeff.get_coeff(A=A1, b=b1, E=E1, T=t)
    k3 = coeff.get_coeff(A=A3, E=E3, T=t)
    k = np.array([k1, k2, k3])
    reaction_rate += [dict(zip(reaction_component, rate_calculator.reaction_rate(v1, v2, x, k)))]

for i in range(len(T)):
    print("under the temperature {}".format(T[i]))
    for (k, v) in reaction_rate[i].items():
        print("the element ", k, "has the reaction rate ", v)