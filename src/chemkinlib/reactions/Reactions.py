
"""Module for setting up different types of reactions."""

import numpy

from chemkinlib.reactions import ReactionRateCoeffs


class ReactionError(Exception):
    """Class for Reaction-related errors."""
    pass

class Reaction(object):
    """Base class for an elementary reaction.
    NOTE: This class is meant to serve as a framework for specific types of reactions!"""
    def __init__(self, rxn_type, is_reversible, rxn_equation, species_list,
                 rate_coeffs_components, reactant_stoich_coeffs, product_stoich_coeffs):
        """
        Initializes Reaction.
    
        INPUTS:
        -------
        rxn_type : str
            type of reaction (e.g. "Elementary")
        is_reversible : bool
            True if reaction is reversible
        rate_coeffs_components : dict
            dictionary of components (e.g. 'A', 'b', and/or 'E')
            to compute reaction rate coefficients
        rxn_equation : str
            string representation of reaction equation
        reactant_stoich_coeffs : dict
            dictionary of integers for reactant stoichiometric coefficients
        product_stoich_coeffs : dict
            dictionary of integers for product stoichiometric coefficients
        unique_species : list
            list of unique chemical species from original xml file (could be a bunch
            of reactions sometimes containing common species)
        species_list : list
            list of chemical species from original xml file (useful for ordering)

        ATTRIBUTES:
        -----------
        temperature : int or float
            temperature of reaction, in Kelvin
        concentrations : list
            concentrations of species involved in reaction
        rxn_rate_coeff : float
            reaction rate coefficient
        """
        self.rxn_type = rxn_type
        self.is_reversible = is_reversible
        self.rate_coeffs_components = rate_coeffs_components
        self.rxn_equation = rxn_equation

        self.reactant_stoich_coeffs = reactant_stoich_coeffs
        self.product_stoich_coeffs = product_stoich_coeffs

        self.unique_species = self.get_unique_species()
        self.species_list = species_list

        # Pad the "nonactive" (non-participating) species with coefficient 0
        for specie in self.species_list:
            if specie not in self.reactant_stoich_coeffs:
                self.reactant_stoich_coeffs[specie] = 0
            if specie not in self.product_stoich_coeffs:
                self.product_stoich_coeffs[specie] = 0
         
        self.temperature = None
        self.concentrations = []
        self.rxn_rate_coeff = None

    def __str__(self):
        """Returns user-friendly string representation of reaction.

        RETURNS:
        --------
        info : str
            string representation of reaction (reaction equation)
        """
        info = "Reaction : {}".format(self.rxn_equation)
        return info

    def __len__(self):
        """Returns number of unique species in reaction.

        RETURNS:
        --------
        n_species : int
            Number of unique species involved in the reaction
        """
        n_species = len(self.unique_species)
        return n_species

    def get_unique_species(self):
        """Helper function to return unique species involved
        in the reaction.
        
        RETURNS:
        --------
        unique_species : list
            list of unique species in reaction
        """
        reactant_species = self.reactant_stoich_coeffs.keys()
        product_species = self.product_stoich_coeffs.keys()
        unique_species = list(set(reactant_species) | set(product_species))
        return unique_species

    def set_temperature(self, T):
        """Sets temperature of the reaction

        INPUTS:
        -------
        T : float
            Temperature of reaction

        NOTES:
        ------
        POST:
            - Updates self.temperature
            - Raises ValueError if inputed temperature is non-positive
        """
        if T <= 0:
            raise ValueError("Temperature has to be a positive value!")

        self.temperature = T

    def set_concentrations(self, X):
        """Sets concentrations of the reaction from a dictionary

        INPUTS:
        -------
        X : dict 
            dictionary with species and corresponding concentrations

        NOTES:
        ------
        PRE:
            - Raises KeyError if input dictionary of species and concentrations
                contains name of species that was not in original xml file
            - Raises ValueError if any of concentrations are negative
        """
        try:
            ordered_concentrations = self.order_dictionaries(X)
        except KeyError:
            raise KeyError("Invalid concentration entered!")

        if (numpy.array(ordered_concentrations) < 0).any():
            raise ValueError("You cannot have negative concentrations!")
        
        self.concentrations = numpy.array(ordered_concentrations)

    def set_concentrations_from_array(self, c):
        """Set the concentrations of the reaction from an array

        :param c:
        :return: None
        """
        self.concentrations = c

    def order_dictionaries(self, dictionary):
        """Helper function to order dictionaries (of concentrations,
        stoichiometric coefficients) based on ordering from species_list.
        This is to ensure a consistent ordering scheme.

        INPUTS:
        -------
        dictionary : dict
            dictionary to order 

        RETURNS:
        --------
        list_of_interest : list
            list of dictionary's keys in order of species_list
        """
        index_map = {v: i for i, v in enumerate(self.species_list)}
        sorted_tuple_list = sorted(dictionary.items(), key=lambda pair: index_map[pair[0]])
        list_of_interest = [element[1] for element in sorted_tuple_list]
        return list_of_interest

    def compute_reaction_rate_coeff(self, T=None):
        """Computes reaction rate coefficients of reaction.

        RETURNS:
        --------
        k : numeric type (or list of numeric type)
            Reaction rate coefficient

        NOTES:
        ------
        POST:
            - Raises NotImplementedError (user must define this function)
        """
        raise NotImplementedError

    def compute_progress_rate(self, T=None):
        """Computes progress rates of reaction.

        RETURNS:
        --------
        omega_array : numpy.ndarray
            Array of progress rates of reaction

        NOTES:
        ------
        POST:
            - Raises NotImplementedError (user must define this function)
        """
        raise NotImplementedError
        
    def compute_reaction_rate(self, T=None):
        """Computes reaction rate of this SINGLE reaction object.

        RETURNS:
        --------
        rxn_rate_array: numpy.ndarray
            Array of reaction rates of reaction

        NOTES
        -----
        POST:
            - Raises ValueError if reactant or product stoichiometric coefficients
                are negative
            - Raises ReactionError if compute_reaction_rate_coeff nor compute_progress_rate
                has not been implemented. This Reaction class is meant to serve as a framework
                for specific types of reactions!
        """
        reactant_stoich_coeffs = numpy.array(self.order_dictionaries(self.reactant_stoich_coeffs))
        product_stoich_coeffs = numpy.array(self.order_dictionaries(self.product_stoich_coeffs))
        concen_array = self.concentrations

        if (reactant_stoich_coeffs < 0).any():
            raise ValueError("Reactant stoichiometric coefficients must be positive!")
        
        if (product_stoich_coeffs < 0).any():
            raise ValueError("Product stoichiometric coefficients must be positive!")

        try:
            progress_rate = self.compute_progress_rate(T)
            nu_i = product_stoich_coeffs - reactant_stoich_coeffs

            reaction_rate_1_eq = progress_rate * nu_i
            return reaction_rate_1_eq

        except NotImplementedError:
            raise ReactionError('''You must first implement the functions to
                                compute the reaction rate coefficients and progress rates!''')




class IrreversibleReactionError(Exception):
    """Error for misclassified IrreversibleReaction."""
    pass

class IrreversibleReaction(Reaction):
    """Class for irreversible elementary reaction"""
    def __init__(self, rxn_type, is_reversible, rxn_equation, species_list, rate_coeffs_components,
                 reactant_stoich_coeffs, product_stoich_coeffs):
        """Initializes IrreversibleReaction.

        NOTES:
        ------
        PRE:
            - Raises IrreversibleReactionError if reaction type is not irreversible OR if reaction
                is not elementary (must satisfy both!)
        """
        super(IrreversibleReaction, self).__init__(rxn_type, is_reversible, rxn_equation,
                                                   species_list, rate_coeffs_components,
                                                   reactant_stoich_coeffs, product_stoich_coeffs)

        if not (rxn_type == "Elementary" and is_reversible == False):
            raise IrreversibleReactionError("This reaction is not irreversible nor elementary!") 

    def compute_reaction_rate_coeff(self, T=None):
        """Computes reaction rate coefficients of reaction.

        INPUTS:
        -------
        T : float
            temperature of reaction, in K

        RETURNS:
        --------
        k : numeric type (or list of numeric type)
            Reaction rate coefficient
        """
        T = self.temperature
        k = ReactionRateCoeffs.ReactionCoeff(self.rate_coeffs_components,
                          T=self.temperature).k
        self.rxn_rate_coeff = k
        return k

    def compute_progress_rate(self, T=None):
        """Computes progress rates of reaction.

        INPUTS:
        -------
        T : float
            temperature of reaction, in K

        RETURNS:
        --------
        omega_array : numpy.ndarray
            Array of progress rates of reaction

        NOTES:
        ------
        PRE:
            - Raises ValueError if concentrations have not been set
        """
        T = self.temperature
        reactant_stoich_coeffs = numpy.array(self.order_dictionaries(self.reactant_stoich_coeffs))
        concen_array = self.concentrations

        if len(concen_array) == 0:
            raise ValueError("You must set the concentrations first!")

        k = self.compute_reaction_rate_coeff(T)

        concen_powered_j = concen_array ** reactant_stoich_coeffs

        progress_rate = k * numpy.prod(concen_powered_j)
        return progress_rate


class ReversibleReactionError(Exception):
    """Error for misclassified ReversibleReaction."""
    pass

class ReversibleReaction(Reaction):
    """Class for reversible reaction"""
    def __init__(self, rxn_type, is_reversible, rxn_equation, species_list, rate_coeffs_components,
                 reactant_stoich_coeffs, product_stoich_coeffs):
        """Initializes ReversibleReaction.

        NOTES:
        ------
        PRE:
            - Raises ReversibleReactionError if reaction type is not reversible OR if reaction
                is not elementary (must satisfy both!)
        """
        super(ReversibleReaction, self).__init__(rxn_type, is_reversible, rxn_equation, species_list, rate_coeffs_components,
                 reactant_stoich_coeffs, product_stoich_coeffs)

        self.NASA_poly_coefs_dict = None
        self.NASA_poly_coefs = None

        if not (rxn_type == "Elementary" and is_reversible == True):
            raise ReversibleReactionError("This reaction is not reversible nor elementary!") 

    def compute_reaction_rate_coeff(self, T=None):
        """Computes reaction rate coefficients of reaction.

        INPUTS:
        -------
        T : float
            temperature of reaction, in K

        RETURNS:
        --------
        kf : numeric type (or list of numeric type)
            forward reaction rate coefficient
        kb : numeric type (or list of numeric type)
            backward reaction rate coefficient

        NOTES:
        ------
        PRE:
            - Raises ValueError if NASA polynomial coefficients have not been
                set before trying to compute reaction rate coefficients
                (See class function set_NASA_poly_coefs())
        """
        coeffs = ReactionRateCoeffs.ReactionCoeff(self.rate_coeffs_components, T=self.temperature)
        self.forward_rxn_rate_coeff = coeffs.k
        
        reactant_stoich_coeffs = numpy.array(self.order_dictionaries(self.reactant_stoich_coeffs))
        product_stoich_coeffs = numpy.array(self.order_dictionaries(self.product_stoich_coeffs))
        nui = product_stoich_coeffs - reactant_stoich_coeffs
        
        if (self.NASA_poly_coefs is None):
            raise ValueError("Must set NASA polynomial coefficients before computing rxn rate coefficients!")
        
        back_coeffs = ReactionRateCoeffs.BackwardCoeff(nui, self.NASA_poly_coefs)
        self.backward_rxn_rate_coeff = back_coeffs.compute_backward_coeffs(self.forward_rxn_rate_coeff,
                                                                           self.temperature)
        return self.forward_rxn_rate_coeff, self.backward_rxn_rate_coeff

    def set_NASA_poly_coefs(self, coefs):
        """Sets NASA polynomial coefficients.
        
        INPUTS:
        -------
        coefs : dict
            dictionary of NASA polynomial coefficients
            (For format, see Parser class' get_NASA_poly_coefs() function)
        """
        # order NASA polynomial coefficients by species (order specified by species_list)
        index_map = {v: i for i, v in enumerate(self.species_list)}
        sorted_tuple_list = sorted(coefs.items(), key=lambda pair: index_map[pair[0]])
        list_of_interest = [element[1] for element in sorted_tuple_list]
        
        # lists of high and low T polynomial coefficients IN ORDER of species
        list_of_interest = numpy.array(list_of_interest)

        self.NASA_poly_coefs_dict = coefs
        self.NASA_poly_coefs = list_of_interest

    def compute_progress_rate(self, T=None):
        """Computes progress rates of reaction.

        INPUTS:
        -------
        T : float
            temperature of reaction, in K

        RETURNS:
        --------
        omega_array : numpy.ndarray
            Array of progress rates of reaction

        NOTES:
        ------
        PRE:
            - Raises ValueError if concentrations have not been not set
        """

        reactant_stoich_coeffs = numpy.array(self.order_dictionaries(self.reactant_stoich_coeffs))
        product_stoich_coeffs = numpy.array(self.order_dictionaries(self.product_stoich_coeffs))
        concen_array = self.concentrations

        if len(concen_array) == 0:
            raise ValueError("You must set the concentrations first!")

        #compute the forward and backward reaction coeffs
        self.compute_reaction_rate_coeff(T)
  
        concen_powered_j_forward = concen_array ** reactant_stoich_coeffs
        concen_powered_j_backward = concen_array ** product_stoich_coeffs

        #compute the forward progress rate
        progress_rate_forward = self.forward_rxn_rate_coeff * numpy.prod(concen_powered_j_forward)

        #compute the backward progress rate
        progress_rate_backward = self.backward_rxn_rate_coeff * numpy.prod(concen_powered_j_backward)

        #code for debugging, save for future use
        #print("progress_rate_forward", progress_rate_forward, "progress_rate_backward", progress_rate_backward)
        #print("concen_powered_j_forward", concen_powered_j_forward, "concen_powered_j_backward", concen_powered_j_backward)
        #print("concen_array", concen_array, "reactant_stoich_coeffs", reactant_stoich_coeffs, "product_stoich_coeffs", product_stoich_coeffs)

        return progress_rate_forward - progress_rate_backward
