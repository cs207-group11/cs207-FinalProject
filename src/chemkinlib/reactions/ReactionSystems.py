#!/usr/bin/env python
# -*- coding: utf-8 -*- 

"""Module for setting up reaction system."""

import numpy

from chemkinlib.reactions import Reactions, ReactionRateCoeffs

from scipy.integrate import ode

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
        self.vis_concentrations = concentrations

        # Set up each reaction
        for r in self.reaction_list:
            r.set_concentrations(concentrations)
            r.set_temperature(self.temperature)
            
            if isinstance(r, Reactions.ReversibleReaction):
                r.set_NASA_poly_coefs(self.NASA_matrix)


        self.concentrations = reaction_list[0].concentrations
        # ODE integrator possible choices of solver: dopri5, dop853(both are explicit ruggi-kutta method) vode, zvode
        self.r = ode(self.compute_reaction_rate).set_integrator("dop853")
        self.r.set_initial_value(self.concentrations, 0)


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
        """
        sort the reaction rate and put it into dictionary format

        RETURNS:
        --------
        rxn_rates_dict: The sorted dictionary of reaction rate
        """
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


    def compute_reaction_rate(self, a, concentrations):
        '''
        Ordinary differential equation for the ODE solver. It gets the current concentration,
        modifies the concentration in each reaction and computes the current reaction rate(first-order derivative)
        INPUTS:
        -------
        a : float
            placeholder of time to match the format of the ode solver.
            In our ODE the derivative doesn't depend on the current time, but only on the current state(concentration)
        concentration: list[float]
            current state (current concentration)

        RETURNS:
        --------
        self.get_reaction_rate() : tuple
            current reaction rate
        '''
        #update the concentration of the species in each reaction and also in the system
        for r in self.reaction_list:
            r.set_concentrations_from_array(concentrations)
        self.concentrations = concentrations
        return self.get_reaction_rate()


    def step(self, dt):
        """Solve the ODE， get the state after dt time

        INPUTS:
        -------
        dt : float
            timestep of the next state

        RETURNS:
        --------
        (self.r.t, self.r.y) : tuple
            current time and current concentration"""
        self.r.integrate(self.r.t + dt)
        return self.r.t, self.r.y
