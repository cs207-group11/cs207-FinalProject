
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
        """
        self.rxn_type = rxn_type
        self.is_reversible = is_reversible
        self.rate_coeffs_components = rate_coeffs_components
        self.rxn_equation = rxn_equation
        self.reactant_stoich_coeffs = reactant_stoich_coeffs
        self.product_stoich_coeffs = product_stoich_coeffs

        self.temperature = None
        self.concentrations = None

    def __str__(self):
        """Returns string representation of reaction."""
        #return self.rxn_equation
        pass # TODO

    def setTemperature(self, T):
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

    def setConcentrations(self, X):
        """Sets concentrations of the reaction

        INPUTS
        ======
        X : 
        """
        # X >= 0, correct dimensions
        self.concentrations = X

    def compute_reaction_rates(self):
        # rxn_rate_coeff = ReactionCoeff()
        # return rxn_rate_coeff.get_coeff(self.rate_coeffs_components,
        #                                 self.temperature)
        pass












