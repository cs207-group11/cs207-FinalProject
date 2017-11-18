
"""Module for parser for reactions."""

from chemkin import *

import sqlite3
import numpy
import os.path
import xml.etree.ElementTree as ET


class ReactionParser():
    """Class for parsing input xml file describing reactions"""

    def __init__(self, xml_filename):
        """Initializes ReactionParser

        INPUTS:
        -------
        xml_filename : str
            filename of input xml file

        ATTRIBUTES:
        -----------
        reaction_list : list
            list of Reaction (or Reaction-inherited) objects

        NOTES:
        ------
        POST:
            - Raises IOError if inputed xml file not found

        """
        if os.path.isfile(xml_filename):
            self.xml_filename = xml_filename
            tree = ET.parse(self.xml_filename)
            self.rxns = tree.getroot()
        else:
            raise IOError("Reaction (xml) file not found!")

        self.reaction_list = []
        self.get_species()
        self.get_NASA_poly_coefs()
        self.get_reaction_list()

    def get_species(self):
        """Returns reaction species
        
        RETURNS:
        --------
        species: a list of species as dictionary
                in the form: key = species name, value = None,
                value will be filled after we obtain the concentration from user.
        """

        phase = self.rxns.findall('phase')
        species_list = []
        for phase in self.rxns.findall('phase'):
            species = phase.find('speciesArray').text
        species_list = species.split()
        self.species = {}
        for specie in species_list:
            self.species[specie] = None
        return self.species

    def get_rxn_type(self, reaction):

        """Returns reaction type

        INPUTS:
        -------
        reaction: parsed xml file which contains information about reactions


        RETURNS:
        --------

        rxn_type: a string describe reaction's type, such as "elementary".
        """
        rxn_type = reaction.get('type')
        return rxn_type

    def get_is_reversible(self, reaction):

        """Returns information about whether the reaction is reversible
        
        INPUTS:
        -------
        reaction: parsed xml file which contains information about reactions
        
        RETURNS:
        --------
        is_reversible: a boolean, True = reversible and False = irreversible
        """

        if reaction.get('reversible') == "yes":
            is_reversible = True
        else:
            is_reversible = False
        return is_reversible


    def get_rxn_equation(self,reaction):
        """Returns reaction equation
        
        INPUTS:
        -------
        reaction: parsed xml file which contains information about reactions
        
        RETURNS:
        --------
        rxn_equation: a string, reaction equation
        """
        rxn_equation = reaction.find('equation').text
        return rxn_equation

    def get_rate_coeffs_components(self,reaction):
        """Returns reaction rate coefficient components
        based on type of coefficient
        
        INPUTS:
        -------
        reaction: parsed xml file which contains information about reactions
        
        RETURNS:
        --------
        rate_coeffs_components: a dictionary, as the form of {coefficient name: coefficient value}. 
        """
        for coef in reaction.findall('rateCoeff'):

            # Arrhenius-type

            if coef.find('Arrhenius') is not None:
                for arr in coef.findall('Arrhenius'):
                    if arr.find('A') is None or arr.find('A').text is None:
                        raise ValueError("Didn't provide coefficient A for modified Arrhenius ")
                    else:
                        A = float(arr.find('A').text)
                    if arr.find('E') is None or arr.find('E').text is None:
                        raise ValueError("Didn't provide coefficient E for modified Arrhenius ")
                    else:
                        E = float(arr.find('E').text)
                rate_coeffs_components = {
                    "A":A,
                    "E":E
                }

            # modified Arrhenius-type

            if coef.find('modifiedArrhenius') is not None:
                for arr in coef.findall('modifiedArrhenius'):
                    if arr.find('A') is None or arr.find('A').text is None:
                        raise ValueError("Didn't provide coefficient A for modified Arrhenius ")
                    else:
                        A = float(arr.find('A').text)
                    if arr.find('b') is None or arr.find('b').text is None:
                        raise ValueError("Didn't provide coefficient b for modified Arrhenius ")
                    else:
                        b = float(arr.find('b').text)
                    if arr.find('E') is None or arr.find('E').text is None:
                        raise ValueError("Didn't provide coefficient E for modified Arrhenius ")
                    else:
                        E = float(arr.find('E').text)
                rate_coeffs_components = {

                    "A":A,
                    "b":b,
                    "E":E
                }

            # constant-type
            if coef.find('Constant') is not None:
                for arr in coef.findall('Constant'):
                    if arr.find('k') is None or arr.find('k').text is None:
                        raise ValueError("Didn't provide coefficient k for modified Arrhenius ")
                    else:
                        k = float(arr.find('k').text)
                rate_coeffs_components = {
                    "k":k
                }
            return rate_coeffs_components

    def get_reactant_stoich_coeffs(self,reaction):
        """Returns reactant stoichiometric coefficients
        
        INPUTS:
        -------
        reaction: parsed xml file which contains information about reactions
        
        RETURNS:
        --------
        reactant_stoich_coeffs: a dictionary, as the form of {reactant name: coefficient}. 
        """
        reactant_stoich_coeffs = {}
        for reactant in reaction.find('reactants').text.split():
            key = reactant.split(":")[0]
            value = reactant.split(":")[1]
            reactant_stoich_coeffs[key] = int(value)
        return reactant_stoich_coeffs
    
    def get_product_stoich_coeffs(self, reaction):
        """Returns product stoichiometric coefficients
        
        INPUTS:
        -------
        reaction: parsed xml file which contains information about reactions
        
        RETURNS:
        --------
        product_stoich_coeffs: a dictionary, as the form of {product name: coefficient}. 
        """
        product_stoich_coeffs = {}
        for product in reaction.find('products').text.split():
            key = product.split(":")[0]
            value = product.split(":")[1]
            product_stoich_coeffs[key] = int(value)
        return product_stoich_coeffs

    def get_coeffs(self, species_name, temp_range):
        """Returns NASA polynomial coefficients from sql database
        from appropriate temperature range.

        INPUTS:
        -------
        species_name : str
            name of chemical species
        temp_range : str
            temperature range for reaction (options: 'low' or 'high')

        RETURNS:
        --------
        nasa_coeffs : numpy.ndarray
            nasa coefficients for species in reaction
        """
        db = sqlite3.connect('NASA_poly_coeff.sqlite')
        cursor = db.cursor()

        if temp_range not in ['high', 'low']:
            raise ValueError("Temperature range can only be 'high' or 'low'...")

        if temp_range == "high":
            cursor.execute('''SELECT COEFF_1, COEFF_2, COEFF_3, COEFF_4, 
                           COEFF_5, COEFF_6, COEFF_7 FROM HIGH WHERE 
                           SPECIES_NAME = ?''', (species_name,))
        
        else:
            cursor.execute('''SELECT COEFF_1, COEFF_2, COEFF_3, COEFF_4, 
                           COEFF_5,COEFF_6,COEFF_7 FROM LOW 
                           WHERE SPECIES_NAME = ?''', (species_name,))
        return numpy.array(cursor.fetchone())

    def get_Tmid(self, species_name):
        """Returns middle T value.

        INPUTS:
        -------
        species_name : str
            name of chemical species

        RETURNS:
        --------
        T_mid : float
            mid T value
        """
        db = sqlite3.connect('NASA_poly_coeff.sqlite')
        cursor = db.cursor()
        cursor.execute('''SELECT TLOW FROM HIGH WHERE SPECIES_NAME = ?''', (species_name,))
        return cursor.fetchone()[0]

    def get_NASA_poly_coefs(self):
        """Parses all information for all reactions

        RETURNS:
        --------
        NASA_poly_coeffs: list[dict]
            List of stored NASA coefficients for each species. Within each dictionary,
                key = low, value = coeff in low temp range and
                key = high, value = coeff in high temp range
        """
        NASA_poly_coefs = []
        for species_name in self.species:
            coef = {}
            coef['Tmid'] = self.get_Tmid(species_name)
            coef['low'] = self.get_coeffs(species_name, 'low')
            coef['high'] = self.get_coeffs(species_name, 'high')
            NASA_poly_coefs.append(coef)
        self.NASA_poly_coefs = NASA_poly_coefs
        return self.NASA_poly_coefs

    def get_reaction_list(self):
        """Appends a Reaction/Reaction-inherited object fo each reaction described by
        input xml file to self.reaction_list.
        """
        for reactionData in self.rxns.findall('reactionData'):
            for reaction in reactionData.findall('reaction'):
                species = self.get_species()
                is_reversible = self.get_is_reversible(reaction) 
                rxn_type = self.get_rxn_type(reaction)
                rxn_equation = self.get_rxn_equation(reaction)
                rate_coeffs_components = self.get_rate_coeffs_components(reaction)
                reactant_stoich_coeffs = self.get_reactant_stoich_coeffs(reaction)
                product_stoich_coeffs = self.get_product_stoich_coeffs(reaction)
                
                # IRREVERSIBLE elementary reaction case
                if is_reversible == False and rxn_type == "Elementary":
                    rxn = IrreversibleReaction(rxn_type, is_reversible, rxn_equation,
                                           species, rate_coeffs_components,
                                           reactant_stoich_coeffs, product_stoich_coeffs)
                    self.reaction_list.append(rxn)

                # REVERSIBLE elementary reaction case
                elif is_reversible == True and rxn_type == "Elementary":
                    rxn = ReversibleReaction(rxn_type, is_reversible, rxn_equation,
                                         species, rate_coeffs_components,
                                         reactant_stoich_coeffs, product_stoich_coeffs)
                    self.reaction_list.append(rxn)

                # Unhandled reaction case
                else:
                    raise NotImplementedError("This type of reaction has not been implemented yet!")
