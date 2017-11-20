
"""Module with classes for computing reaction rate coefficients (forward and backward)."""

import numbers
import numpy
import warnings

class ReactionCoeff(object):
    """Class for reaction rate coefficients, or values k."""
    def __init__(self, k_parameters, T=None):
        """Initializes reaction rate coefficients.

        INPUTS:
        -------
        T : int or float
            temperature of the reaction (in Kelvin)
        k_parameters : dictionary
            dictionary of parameters to compute k
        """
        self.k_parameters = k_parameters
        self.T = T
        self.k = self.get_coeff(self.k_parameters, self.T)

    def get_coeff(self, k_parameters, T):
        """Computes reaction rate coefficients depending on passed parameters.

        INPUTS:
        -------
        T : int or float
            temperature of the reaction (in Kelvin)
        k_parameters : dictionary
            dictionary of parameters to compute k
        
        RETURNS:
        --------
        k : int or float
            reaction rate coefficient of the reaction

        NOTES:
        ------
        PRE:
            - Raise ValueError if customized reaction rate coefficient depends on T
        POST:
            - Raises NotImplementedError if dictionary of k parameters is not recognized
            - Options to alter values of R (to change units) but strongly discouraged
            - Raises ValueError if valid T not inputed/set for Arrhenius and modified Arrhenius
        """
        #check if the key-arguments are valid
        keys = set(k_parameters.keys())
        valid_keys = set(['A', 'E', 'b', "R", "k"])
        if not(keys <= valid_keys):
            raise ValueError("Invalid key in the input. Use get_coeff function to implement your own k!")

        # Constant
        if "k" in k_parameters:
            return self.const(k_parameters['k'])
        
        # Arrhenius
        elif ("A" in k_parameters and "E" in k_parameters and
              "b" not in k_parameters):

            if T == None:
                raise ValueError("Temperature has not been set in the reaction. Please use set_temperature() method.")

            if "R" in k_parameters:
                return self.arr(A=k_parameters['A'],
                                E=k_parameters['E'],
                                T=T,
                                R=k_parameters['R'])
            else:
                return self.arr(A=k_parameters['A'],
                                E=k_parameters['E'],
                                T=T)
        
        # Modified Arrhenius
        elif ("A" in k_parameters and "E" in k_parameters  and "b" in k_parameters):
            if T == None:
                raise ValueError("Temperature has not been set in the reaction. Please use set_temperature() method.")

            if "R" in k_parameters:
                return self.mod_arr(A=k_parameters['A'],
                                    E=k_parameters['E'],
                                    R=k_parameters['R'],
                                    b=k_parameters['b'],
                                    T=T)
            else:
                return self.mod_arr(A=k_parameters['A'],
                                    E=k_parameters['E'],
                                    b=k_parameters['b'],
                                    T=T)

        else:
            raise NotImplementedError("The combination of parameters entered is not supported for the calculation of Reaction Rate Coefficient.")


    def const(self, k):
        """Returns constant reaction rate coefficients k.

        INPUTS:
        -------
        k : numeric type 
            constant reaction rate coefficient

        RETURNS:
        --------
        k : numeric type
            constant reaction rate coefficients.

        NOTES:
        ------
        POST:
            - Raises ValueError if k is non-positive!
        """
        if k <= 0:
            raise ValueError("Reaction rate must be positive.")
        
        return k

    def arr(self, A, E, T, R=8.314):
        """Returns Arrhenius reaction rate coefficients k.

        INPUTS:
        -------
        A: float, strictly positive, no default value
           The Arrhenius prefactor
        E: float, no default value
           The Arrhenius parameter
        T: float, strictly positive, no default value
           Temperature T, asuuming a Kelvin scale
        R: float, default value is 8.314, cannot be changed except to convert units
           The ideal gas constant

        RETURNS:
        --------
        k: Arrhenius reaction rate coefficients k,
           floats
           unless A or T is not postive
           in which case a ValueError exception is raised

        NOTES:
        ------
        POST:
            - Raises ValueError if A, T, or R is non-positive
            - Raises Warning if user changes value of R
        """
        if (A <= 0):
            raise ValueError("Arrhenius prefactor 'A' must be positive.")

        if (T <= 0):
            raise ValueError("Temperatures 'T' must be positive.")

        if (R <= 0):
            raise ValueError("Gas constant 'R' must be positive.")

        if not numpy.isclose(R, 8.314):
            warnings.warn("Please do not change the value of"
                          " Universal Gas Constant 'R' unless you are converting units.")

        k = A * numpy.exp(-E / R / T)
        return k

    def mod_arr(self, A, b, E, T, R=8.314):
        """Returns Arrhenius reaction rate coefficients k.

        INPUTS:
        -------
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

        RETURNS:
        --------
        k: Arrhenius reaction rate coefficients k,
           floats
           unless A or T is not postive
           in which case a ValueError exception is raised
           Or b is not a real number
           in which case a TypeError exception is raised

        NOTES:
        ------
        POST:
            - Raises ValueError if A, T, or R is non-positive
            - Raises TypeError if b is not real
            - Raises Warning if user changes value of R
        """

        if (A <= 0):
            raise ValueError("Parameter 'A' must be positive.")
        
        if (T <= 0):
            raise ValueError("Parameter 'T' must be positive.")
        
        if (isinstance(b, numbers.Real)) == False:
            raise TypeError("Parameter 'b' must be a real number.")
        
        if (R <= 0):
            raise ValueError("Gas constant 'R' must be positive.")

        if not numpy.isclose(R, 8.314):
            warnings.warn("Please do not change the value of"
                          " Universal Gas Constant 'R' unless for converting units.")

        k = A * T ** b * numpy.exp(-E / R / T)
        return k




class BackwardCoeff():
    """Class for computing backward reaction rate
    coefficients for reversible reactions."""
    def __init__(self, nui, nasa7_coeffs):
        """Initializes BackwardCoeff.

        INPUTS:
        -------
        nui : numpy.ndarray
            stoichiometric coefficient difference (stoich_products - stoich_reactants)
                for a single reversible reaction
        nasa7_coeffs : numpy.ndarray
            NASA polynomial coefficients (from appropriate temperature range)
                corresponding to species in reversible reaction

        ATTRIBUTES:
        -----------
        p0 : float
            pressure of reaction, in Pascals
        R : float
            gas constant, in J / mol / K
        gamma : numpy.ndarray
            sum of stoichiometric coefficient difference 
        """
        self.nui = nui
        self.nasa7_coeffs = nasa7_coeffs
        self.p0 = 1.0e+05
        self.R = 8.3144598
        self.gamma = numpy.sum(self.nui)

    # def Cp_over_R(self, T):
    #     """Returns specific heat of each specie given by
    #     the NASA polynomials.

    #     INPUTS:
    #     -------
    #     T : float
    #         temperature of reaction

    #     RETURNS:
    #     --------
    #     Cp_R : numpy.ndarray
    #         specific heat values for each specie
    #     """
    #     a = self.nasa7_coeffs
    #     Cp_R = (a[:, 0] + a[:, 1] * T + a[:, 2] * T ** 2.0
    #             + a[:, 3] * T ** 3.0 + a[:, 4] * T ** 4.0)
    #     return Cp_R

    def H_over_RT(self, T):
        """Returns the enthalpy of each specie given by
        the NASA polynomials.

        INPUTS:
        -------
        T : float
            temperature of reaction

        RETURNS:
        --------
        H_RT : numpy.ndarray
            enthalpy values for each specie

        NOTES:
        ------
        PRE:
            - Raises ValueError if inputed temperature is non-positive
        """
        if T <= 0:
            raise ValueError("Temperature has to be a positive value!")

        a = self.nasa7_coeffs
        H_RT = (a[:, 0] + (0.5 * a[:, 1] * T) + (a[:, 2] * T ** 2.0) / 3.0
                + (a[:, 3] * T ** 3.0) / 4.0 + (a[:, 4] * T ** 4.0) / 5.0
                + a[:, 5] / T)
        return H_RT

    def S_over_R(self, T):
        """Returns the entropy of each specie given by
        the NASA polynomials.

        INPUTS:
        -------
        T : float
            temperature of reaction

        RETURNS:
        --------
        S_R : numpy.ndarray
            entropy values for each specie

        NOTES:
        ------
        PRE:
            - Raises ValueError if inputed temperature is non-positive
        """
        if T <= 0:
            raise ValueError("Temperature has to be a positive value!")

        a = self.nasa7_coeffs
        S_R = (a[:, 0] * numpy.log(T) + a[:, 1] * T + (a[:, 2] * T ** 2.0) / 2.0
               + (a[:, 3] * T ** 3.0) / 3.0 + (a[:, 4] * T ** 4.0) / 4.0 + a[:, 6])
        return S_R

    def compute_backward_coeffs(self, kf, T):
        """Returns the backward reaction rate
        coefficient for each specie.

        INPUTS:
        -------
        kf : numpy.ndarray[float]
            array of forward reaction rate coefficients for each specie in reaction
        T : float
            temperature of reaction

        RETURNS:
        --------
        kb : backward reaction rate coefficient for each specie

        NOTES:
        ------
        PRE:
            - Raises ValueError if inputed temperature is non-positive
        """
        if T <= 0:
            raise ValueError("Temperature has to be a positive value!")

        # Change in enthalpy and entropy for each reaction
        delta_H_over_RT = numpy.dot(self.nui, self.H_over_RT(T))
        delta_S_over_R = numpy.dot(self.nui, self.S_over_R(T))

        # Negative of change in Gibbs free energy for each reaction
        delta_G_over_RT = delta_S_over_R - delta_H_over_RT

        # Prefactor in k_e (equilibrium coefficient)
        fact = self.p0 / self.R / T

        ke = (fact ** self.gamma) * numpy.exp(delta_G_over_RT)

        return kf / ke
