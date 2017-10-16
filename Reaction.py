
"""Base class for Reaction."""

import numpy as np


class Reaction():
    """Class for Reaction"""
    def __init__(self, rxn_type, is_reversible,
                 rxn_equation, rate_coeffs_components,
                 reactant_stoich_coeffs, product_stoich_coeffs):
        """Initializes Reaction
    
        INPUTS:
        -------
        rxn_type : str
            type of reaction (e.g. "Elementary")
        is_reversible : bool
            True if reaction is reversible
        rxn_equation : str
            string representation of reaction equation
        rate_coeffs_components : dict
            dictionary of components (e.g. 'A', 'b', and/or 'E')
            to compute reaction rate coefficients
        reactant_stoich_coeffs : dict
            dictionary of integers for reactant stoichiometric coefficients
        product_stoich_coeffs : dict
            dictionary of integers for product stoichiometric coefficients

        ATTRIBUTES:
        -----------
        temperature : int or float
            user-inputed temperature (in K)
        concentrations : list
            list of concentrations of species involved in reaction
        rxn_rate_coeff : float
            reaction rate coefficient (will be computed later)
        """
        self.rxn_type = rxn_type
        self.is_reversible = is_reversible
        self.rate_coeffs_components = rate_coeffs_components
        self.rxn_equation = rxn_equation
        self.reactant_stoich_coeffs = reactant_stoich_coeffs
        self.product_stoich_coeffs = product_stoich_coeffs

        self.temperature = None
        self.concentrations = None
        self.rxn_rate_coeff = None

    def __str__(self):
        """Returns user-friendly string representation of reaction.

        RETURNS
        =======
        info : str
            string representation of reaction (reaction equation)
        """
        info = "Reaction : {}".format(self.rxn_equation)
        return info

    def __len__(self):
        """Returns number of unique species in reaction.

        RETURNS
        =======
        n_species : int
            Number of unique species involved in the reaction
        """
        reactant_species = self.reactant_stoich_coeffs.keys()
        product_species = self.product_stoich_coeffs.keys()
        n_species = len(list(set(reactant_species) | set(product_species)))
        return n_species

    def set_temperature(self, T):
        """Sets temperature of the reaction

        INPUTS 
        ======
        T : float
            Temperature of reaction

        NOTES
        =====
        POST:
            - Updates self.temperature
            - Raises ValueError if inputed temperature is non-positive
        """
        if T <= 0:
            raise ValueError("Temperature has to be a positive value!")
        self.temperature = T

    def set_concentrations(self, X):
        """Sets concentrations of the reaction

        INPUTS
        ======
        X : dict 
            dictionary with species and corresponding concentrations
        """
        # X >= 0, correct dimensions
        self.concentrations = X

    def order_concentrations(self, concentration_dict):
        # take in self.concen and order it (order in species in xml)
        pass

    def compute_reaction_coeff(self):
        rxn_rate_coeff = ReactionCoeff()
        self.rxn_rate_coeff = rxn_rate_coeff.get_coeff(
                                self.rate_coeffs_components,
                                self.temperature)

    # def compute_progress_rate(reactant_stoich_coeffs, concen_array, k):
    #     """Computes progress rate for inputed reaction
    #     INPUTS
    #     ======= 
    #     reactant_stoich_coeffs: numpy.ndarray
    #         Stoichiometric coefficients of reactants
    #     concen_array: numpy.ndarray
    #         Concentrations of species
    #     k : numeric type (or list of numeric type)
    #         Reaction rate coefficient
    #     RETURNS
    #     ========
    #     omega_array: numpy.ndarray
    #         Array of progress rate(s) of reaction
    #     NOTES
    #     ======
    #     PRE:
    #          - reactant_stoich_coeffs, an array of positive ints of shape (i, j)
    #             where i = number of species, j = number of reactions
    #          - concen, an array of positive ints of shape (i, )
    #          - k a numeric type or list (list if k of each reaction differs for multiple elementary reactions)
    #     POST:
    #          - raises ValueErrors if stoichiometric coefficients, concentrations,
    #             or reaction rate constants non-positive
    #          - returns array of floats
    #     EXAMPLES
    #     =========
    #     >>> compute_progress_rate(np.array([2., 1., 0.]), np.array([1., 2., 3.]), 10)
    #     array([ 20.])
    #     >>> compute_progress_rate(np.array([[1.0, 2.0], [2.0, 0.0], [0.0, 2.0]]), np.array([1., 2., 1.]), 10)
    #     array([ 40.,  10.])
    #     >>> compute_progress_rate(np.array([[1.0, 2.0], [2.0, 0.0], [0.0, 2.0]]), np.array([1., 2., 1.]), [50., 30.])
    #     array([ 200.,   30.])
    #     """
    #     if (reactant_stoich_coeffs < 0).any():
    #         raise ValueError("Stoichiometric coefficients must be positive!")

    #     if (concen_array <= 0).any():
    #         raise ValueError("Concentrations must be positive!")

    #     if reactant_stoich_coeffs.shape[0] != len(concen_array):
    #         raise ValueError("Number of species must stay consistent (lengths of concen_array and number of columns in coeff array)")

    #     try:
    #         n_rxns = reactant_stoich_coeffs.shape[1]
    #     except IndexError:
    #         n_rxns = 1

    #     omega_array = np.zeros(n_rxns)

    #     for j in range(n_rxns):
    #         if n_rxns == 1:
    #             concen_powered_j = concen_array**reactant_stoich_coeffs
    #         else:
    #             concen_powered_j = concen_array**reactant_stoich_coeffs[:, j]

    #         if isinstance(k, float) or isinstance(k, int):
    #             if k <= 0:
    #                 raise ValueError("Reaction rate constants must be positive!")
                
    #             omega_j = k * np.prod(concen_powered_j)
    #             omega_array[j] = omega_j

    #         elif isinstance(k, list):
    #             if len(k) != n_rxns:
    #                 raise ValueError("If k is a list, its length must equal the number of elementary reactions!")

    #             if (np.array(k) <= 0).any():
    #                 raise ValueError("Reaction rate constants must be positive!")

    #             omega_j = k[j] * np.prod(concen_powered_j)
    #             omega_array[j] = omega_j

    #     return omega_array


    def compute_reaction_rate(self):

        if (self.rxn_type == "Elementary" and
            not self.is_reversible):

            # omega_array = self.compute_progress_rate(reactant_stoich_coeffs, concen_array, k)
            # nu_ij = product_stoich_coeffs - reactant_stoich_coeffs
            rxn_rate_array = np.dot(nu_ij, omega_array)
            return rxn_rate_array
            
           
            
            pass
            # TODO!

        else:
            raise NotImplementedError("Computing reaction rate for this type and reversibility "
                                      "of reaction has not been considered!")
        return None















