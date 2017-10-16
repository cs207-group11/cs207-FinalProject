
import numpy as np

class ReactionRate():
    """Base class for reaction rates"""
    def __init__(self):
        raise NotImplementedError

    def progress_rate(self):
        raise NotImplementedError

    def reaction_rate(self):
        raise NotImplementedError


class IrreversibleReactionRate(ReactionRate):
    def __init__(self):
        pass

    # def progress_rate(self, v, x, k):
    #     """Returns the progress rate for a system of reactions

    #     INPUTS
    #     =======
    #     v: The matrix of Stoichiometric coefficients of reactants
    #     x: The vector of concentration of species
    #     k: The vector of reaction rate coefficient

    #     RETURNS
    #     ========
    #     w: Progress rate of the reaction system

    #     EXAMPLES
    #     =========
    #     >>> r = IrreversibleReactionRate()
    #     >>> r.progress_rate(np.array([[1,2,0], [2, 0, 2]]), np.array([1,2,1]), np.array([10,10]))
    #     array([ 40.,  10.])
    #     """
    #     w = np.zeros(len(v))
    #     for i in range(len(v)):
    #         w[i] = k[i] * np.prod(np.power(x, v[i]))
    #     return w


    def progress_rate(self, reactant_stoich_coeffs, concen_array, k):
        """Computes progress rate for inputed reaction

        INPUTS
        ======= 
        reactant_stoich_coeffs: numpy.ndarray
            Stoichiometric coefficients of reactants
        concen_array: numpy.ndarray
            Concentrations of species
        k : numeric type (or list of numeric type)
            Reaction rate coefficient

        RETURNS
        ========
        omega_array: numpy.ndarray
            Array of progress rate(s) of reaction

        NOTES
        ======
        PRE:
             - reactant_stoich_coeffs, an array of positive ints of shape (i, j)
                where i = number of species, j = number of reactions
             - concen, an array of positive ints of shape (i, )
             - k a numeric type or list (list if k of each reaction differs for multiple elementary reactions)
        POST:
             - raises ValueErrors if stoichiometric coefficients, concentrations,
                or reaction rate constants non-positive
             - returns array of floats

        EXAMPLES
        =========
        >>> compute_progress_rate(np.array([2., 1., 0.]), np.array([1., 2., 3.]), 10)
        array([ 20.])

        >>> compute_progress_rate(np.array([[1.0, 2.0], [2.0, 0.0], [0.0, 2.0]]), np.array([1., 2., 1.]), 10)
        array([ 40.,  10.])

        >>> compute_progress_rate(np.array([[1.0, 2.0], [2.0, 0.0], [0.0, 2.0]]), np.array([1., 2., 1.]), [50., 30.])
        array([ 200.,   30.])
        """
        if (reactant_stoich_coeffs < 0).any():
            raise ValueError("Stoichiometric coefficients must be positive!")

        if (concen_array <= 0).any():
            raise ValueError("Concentrations must be positive!")

        if reactant_stoich_coeffs.shape[0] != len(concen_array):
            raise ValueError("Number of species must stay consistent (lengths of concen_array and number of columns in coeff array)")

        try:
            n_rxns = reactant_stoich_coeffs.shape[1]
        except IndexError:
            n_rxns = 1

        omega_array = np.zeros(n_rxns)

        for j in range(n_rxns):
            if n_rxns == 1:
                concen_powered_j = concen_array**reactant_stoich_coeffs
            else:
                concen_powered_j = concen_array**reactant_stoich_coeffs[:, j]

            if isinstance(k, float) or isinstance(k, int):
                if k <= 0:
                    raise ValueError("Reaction rate constants must be positive!")
                
                omega_j = k * np.prod(concen_powered_j)
                omega_array[j] = omega_j

            elif isinstance(k, list):
                if len(k) != n_rxns:
                    raise ValueError("If k is a list, its length must equal the number of elementary reactions!")

                if (np.array(k) <= 0).any():
                    raise ValueError("Reaction rate constants must be positive!")

                omega_j = k[j] * np.prod(concen_powered_j)
                omega_array[j] = omega_j

        return omega_array


    def reaction_rate(self, reactant_stoich_coeffs, product_stoich_coeffs, concen_array, k):
        """Computes reaction rate for inputed reaction 

        INPUTS
        ======= 
        reactant_stoich_coeffs: numpy.ndarray
            Stoichiometric coefficients of reactants
        product_stoich_coeffs: numpy.ndarray
            Stoichiometric coefficients of products
        concen_array: numpy.ndarray
            Concentrations of species
        k : numeric type (or list of numeric type)
            Reaction rate coefficient(s)

        RETURNS
        ========
        rxn_rate_array: numpy.ndarray
            Array of reaction rates of reaction

        NOTES
        ======
        PRE:
             - reactant_stoich_coeffs, an array of positive ints of shape (i, j)
                where i = number of species, j = number of reactions
             - product_stoich_coeffs, an array of positive ints of shape (i, j)
             - concen, an array of positive ints of shape (i, )
             - k a numeric type or list (list if k of each reaction differs for multiple elementary reactions)
        POST:
             - raises ValueErrors if stoichiometric coefficients, concentrations,
                or reaction rate constants non-positive
             - raises ValueError if stoichiometric coefficient arrays have different dimensions
             - returns array of floats
        
        EXAMPLES
        =========
        >>> compute_rxn_rate(np.array([[1.0, 0.0], [2.0, 0.0], [0.0, 2.0]]), np.array([[0.0, 1.0], [0.0, 2.0], [1.0, 0.0]]), np.array([1., 2., 1.]), 10)
        array([-30., -60.,  20.])
        
        >>> compute_rxn_rate(np.array([[1.0, 0.0], [2.0, 0.0], [0.0, 2.0]]), np.array([[0.0, 1.0], [0.0, 2.0], [1.0, 0.0]]), np.array([1., 2., 1.]), [20, 10])
        array([ -70., -140.,   60.])
        """
        if ((reactant_stoich_coeffs < 0).any() or
            (product_stoich_coeffs < 0).any()):
            raise ValueError("Stoichiometric coefficients must be positive!")

        if (reactant_stoich_coeffs.shape != product_stoich_coeffs.shape):
            raise ValueError("Dimension mismatch for reactant and product stoichiometric coefficients!")

        if (concen_array <= 0).any():
            raise ValueError("Concentrations must be positive!")

        omega_array = self.progress_rate(reactant_stoich_coeffs, concen_array, k)

### HARDCODED
        omega_array = np.ones(3) * omega_array[0]

        nu_ij = product_stoich_coeffs - reactant_stoich_coeffs

        rxn_rate_array = np.dot(nu_ij, omega_array)
        return rxn_rate_array


class ReversibleReactionRate(ReactionRate):
    def progress_rate(self):
        raise NotImplementedError

    def reaction_rate(self):
        raise NotImplementedError






class ReactionCoeff():
    """Class for reaction rate coefficients."""




    """This module contains functions that return reaction rate coefficients

       MEMBERS:
       ========
       const:   Returns a constant reaction rate coefficient
       arr:     Returns the Arrhenius reaction rate coefficient
       mod_arr: Returns the modified Arrhenius reaction rate coefficient
    """
    def __init__(self):
        pass

    def get_coeff(self, T, k_parameters):
        """
        :param k_parameters: k, A, b, E ,T, R(optional, usually not changed)
        :return: reaction coefficient
        """
        if "k" in k_parameters:
            return self.const(k_parameters['k'])
        if "A" in k_parameters and "E" in k_parameters and "b" not in k_parameters:
            if "R" in k_parameters:
                return self.arr(A = k_parameters['A'], E = k_parameters['E'], T = T, R = k_parameters['R'])
            else:
                return self.arr(A= k_parameters['A'], E=k_parameters['E'], T=T)
        if "A" in k_parameters and "E" in k_parameters  and "b" in k_parameters:
            if "R" in k_parameters:
                return self.mod_arr(A = k_parameters['A'], E = k_parameters['E'], R = k_parameters['R'], b = k_parameters['b'], T=T)
            else:
                return self.mod_arr(A=k_parameters['A'], E=k_parameters['E'], b=k_parameters['b'], T=T)

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


if __name__ == "__main__":


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