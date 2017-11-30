
"""Module for setting up reaction system."""

import numpy

from chemkinlib.reactions import Reactions, ReactionRateCoeffs


class ReactionSystemError(Exception):
    """Class for ReactionSystem-related errors."""
    pass

class ReactionSystem(object):
    """Class for a system of reactions extracted from an xml file"""
    def __init__(self, reaction_list, NASA_poly_coefs, temperature, concentrations):
        """Initializes ReactionSystem.
        
        INPUTS:
        -------
        reaction_list : list[Reaction] or list[Reaction-like]
            list of Reaction or Reaction-inherited objects
        NASA_poly_coefs : 

        temp : numeric type
            temperature of reaction
        concentrations : dict
            dictionary of concentrations with key as species name

        ATTRIBUTES:
        -----------
        involved_species : list[str]
            list of all species in reaction system
        """
        self.reaction_list = reaction_list

        if temperature <= 0:
            raise ValueError("Temperature has to be a positive value!")

        self.temperature = temperature
        self.NASA_matrix = self.get_nasa_matrix(NASA_poly_coefs)

        # NOTE: Not sure if this is good enough!
        self.involved_species = reaction_list[0].species_list

        # Set up each reaction
        for r in self.reaction_list:
            r.set_concentrations(concentrations)
            r.set_temperature(self.temperature)
            
            if isinstance(r, Reactions.ReversibleReaction):
                r.set_NASA_poly_coefs(self.NASA_matrix)
           
    def get_reaction_rate(self):
        """Fetches reaction rate for each reaction.

        RETURNS:
        --------
        list_rxn_rates : list[float]
            list of reaction rates of reactions in the system
        """
        reaction_rate_list = [rxnObj.compute_reaction_rate() for rxnObj in self.reaction_list]
        reaction_rate_list = numpy.array(reaction_rate_list)
        rxnrates = numpy.sum(reaction_rate_list, axis=0)
        return rxnrates

    def sort_reaction_rates(self):
        rxn_rates_dict = {}
        list_species_ordered = list(self.involved_species)
        rxnrate = self.get_reaction_rate()

        for i in range(len(rxnrate)):
            rxn_rates_dict[list_species_ordered[i]] = rxnrate[i]
         
        return rxn_rates_dict

    def get_nasa_matrix(self, NASA_poly_coef):
        """Computes array of NASA polynomial coefficients.

        INPUTS:
        -------
        NASA_poly_coef : list[dict]
            list of dictionaries of NASA polynomial coefficients
                labeled by temperature range

        RETURNS:
        --------
        NASA_array : numpy.ndarray
            array of NASA polynomial coefficients for given temperature range
        """

        NASA = {}
        
        for specie in NASA_poly_coef:
            specie_dict = NASA_poly_coef[specie]
            if self.temperature <= specie_dict["Tmid"]: # get the low temperature
                NASA[specie] = specie_dict["low"]
            else:
                NASA[specie] = specie_dict["high"]
        return NASA
