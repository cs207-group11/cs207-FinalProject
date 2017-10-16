
"""Base class for Reaction."""

import numpy as np
from Scratch_from_Shiyu import *


class Reaction():
    """Base class for a reaction"""
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
        species_list : list
            list of unique chemical species involved in reaction (useful for ordering)
        """
        self.rxn_type = rxn_type
        self.is_reversible = is_reversible
        self.rate_coeffs_components = rate_coeffs_components
        self.rxn_equation = rxn_equation

        self.reactant_stoich_coeffs = reactant_stoich_coeffs
        self.product_stoich_coeffs = product_stoich_coeffs

        self.species_list = self.get_unique_species()

        for specie in self.species_list:
            if specie not in self.reactant_stoich_coeffs:
                self.reactant_stoich_coeffs[specie] = 0
            if specie not in self.product_stoich_coeffs:
                self.product_stoich_coeffs[specie] = 0

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
        n_species = len(self.get_unique_species())
        return n_species

    def get_unique_species(self):
        """Helper function to return unique species involved
        in the reaction.

        RETURNS
        =======
        unique_species : list
            List of unique species in reaction
        """
        reactant_species = self.reactant_stoich_coeffs.keys()
        product_species = self.product_stoich_coeffs.keys()
        unique_species = list(set(reactant_species) | set(product_species))
        return unique_species

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

    def order_dictionaries(self, dictionary):
        """Orders dictionaries (of concentrations, stoichiometric coefficients)
        based on ordering from species list from xml file. This is to 
        ensure a consistent ordering scheme.

        INPUTS
        ======
        dictionary : dict
            dictionary to order 
        """
        index_map = {v: i for i, v in enumerate(self.species_list)}
        sorted_tuple_list = sorted(dictionary.items(), key=lambda pair: index_map[pair[0]])
        list_of_interest = [element[1] for element in sorted_tuple_list]
        return list_of_interest

    def compute_reaction_coeff(self):
        rxn_rate_coeff = ReactionCoeff()
        self.rxn_rate_coeff = rxn_rate_coeff.get_coeff(
                                self.temperature,
                                self.rate_coeffs_components)
        return self.rxn_rate_coeff


    def compute_reaction_rate(self):
        if (self.rxn_type == "Elementary" and
            not self.is_reversible):

            reactant_stoich_coeffs = self.order_dictionaries(self.reactant_stoich_coeffs)
            product_stoich_coeffs = self.order_dictionaries(self.product_stoich_coeffs)
            concen_array = self.order_dictionaries(self.concentrations)
            k = self.compute_reaction_coeff()
            
            reactant_stoich_coeffs = np.array(reactant_stoich_coeffs)
            product_stoich_coeffs = np.array(product_stoich_coeffs)
            concen_array = np.array(concen_array)

            rxn_rate_obj = IrreversibleReactionRate()
            rxn_rate = rxn_rate_obj.reaction_rate(reactant_stoich_coeffs,
                                                  product_stoich_coeffs,
                                                  concen_array, k)


            return rxn_rate

        else:
            raise NotImplementedError("Computing reaction rate for this type and reversibility "
                                      "of reaction has not been considered!")
        return None





if __name__ == "__main__":

    test = Reaction(rxn_type="Elementary",
                    is_reversible=False,
                    rxn_equation="A + B =] C",
                    rate_coeffs_components={'k': 10},
                    reactant_stoich_coeffs={'A' :2, 'B':1, 'C':0},
                    product_stoich_coeffs={'A' :0, 'B':0, 'C':1})

    print(test)



    test.set_temperature(10)
    test.set_concentrations({'A':1, 'B': 2, 'C':3})

    k = test.compute_reaction_coeff()

    rxn_rate_test = test.compute_reaction_rate()
    print(rxn_rate_test)









