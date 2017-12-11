
"""Module for parsing reaction parameters."""

import numpy
import os.path
import sqlite3
import xml.etree.ElementTree as ET

# path to the xml files
from chemkinlib.config import DATA_DIRECTORY
from chemkinlib.reactions import Reactions


class ReactionParser():
    """Class for parsing input xml file describing reactions."""
    def __init__(self, xml_filename):
        """Initializes ReactionParser.

        INPUTS:
        -------
        xml_filename : str
            filename of input xml file

        ATTRIBUTES:
        -----------
        reaction_list : list
            list of Reaction (or Reaction-inherited) objects
        rxns : xml.etree.ElementTree.Element
            root node of xml tree
        species : dict[str]
            dictionary of specie names
        NASA_poly_coefs : dict
            dictionary that will contain NASA polynomial
            coefficients at high and low temperature ranges
            (see class function get_NASA_poly_coefs() for exact
            form)

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
        self.get_species() # updates species
        self.get_NASA_poly_coefs() # updates nasa polynoms
        self.get_reaction_list() # appendsto reaction_list

    def get_species(self):
        """Returns reaction species from input.
        
        RETURNS:
        --------
        species: dict[str]  
            dictionary of the form: key = species name, value = None
            (used to order information by species)
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
        """Returns reaction type from input.

        INPUTS:
        -------
        reaction : xml.etree.ElementTree.Element
            parsed xml file information about reactions

        RETURNS:
        --------
        rxn_type : str
            string describing reaction type (e.g. "elementary")
        """
        rxn_type = reaction.get('type')
        return rxn_type

    def get_is_reversible(self, reaction):
        """Returns information about whether the reaction is reversible.
        
        INPUTS:
        -------
        reaction : xml.etree.ElementTree.Element
            parsed xml file information about reactions
        
        RETURNS:
        --------
        is_reversible : bool
            if True, reversible
            if False, irreversible
        """
        if reaction.get('reversible') == "yes":
            is_reversible = True
        else:
            is_reversible = False
        return is_reversible

    def get_rxn_equation(self,reaction):
        """Returns reaction equation from input.
        
        INPUTS:
        -------
        reaction : xml.etree.ElementTree.Element
            parsed xml file information about reactions
        
        RETURNS:
        --------
        rxn_equation : str
            a string representation of reaction equation
        """
        rxn_equation = reaction.find('equation').text
        return rxn_equation

    def get_rate_coeffs_components(self, reaction):
        """Returns reaction rate coefficient components
        based on type of coefficient.
        
        INPUTS:
        -------
        reaction : xml.etree.ElementTree.Element
            parsed xml file information about reactions
        
        RETURNS:
        --------
        rate_coeffs_components : dict
            a dictionary of the form {coefficient component name: coefficient component value}. 
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

    def get_reactant_stoich_coeffs(self, reaction):
        """Returns reactant stoichiometric coefficients from input.
        
        INPUTS:
        -------
        reaction : xml.etree.ElementTree.Element
            parsed xml file information about reactions
        
        RETURNS:
        --------
        reactant_stoich_coeffs : dict
            dictionary in the form {reactant name: stoich coefficient}. 
        """
        reactant_stoich_coeffs = {}
        for reactant in reaction.find('reactants').text.split():
            key = reactant.split(":")[0]
            value = reactant.split(":")[1]
            reactant_stoich_coeffs[key] = int(value)
        return reactant_stoich_coeffs
    
    def get_product_stoich_coeffs(self, reaction):
        """Returns product stoichiometric coefficients from input.
        
        INPUTS:
        -------
        reaction : xml.etree.ElementTree.Element
            parsed xml file information about reactions
        
        RETURNS:
        --------
        product_stoich_coeffs : dict
            dictionary in the form {product name: stoich coefficient}. 
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
        db = sqlite3.connect(DATA_DIRECTORY + '/' + 'NASA_poly_coeffs.sqlite')
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

        nasa_coeffs = numpy.array(cursor.fetchone())
        db.commit()
        db.close()
        return nasa_coeffs

    def get_Tmid(self, species_name):
        """Returns middle temperature (T) value.

        INPUTS:
        -------
        species_name : str
            name of chemical species

        RETURNS:
        --------
        T_mid : float
            mid temperature value (ie separates
            low and high temperature ranges)
        """
        db = sqlite3.connect(DATA_DIRECTORY + '/' + 'NASA_poly_coeffs.sqlite')
        cursor = db.cursor()
        cursor.execute('''SELECT TLOW FROM HIGH WHERE SPECIES_NAME = ?''', (species_name,))
        Tmid = cursor.fetchone()[0]
        db.commit()
        db.close()
        return Tmid

    def get_NASA_poly_coefs(self):
        """Fetches and updates all NASA polynomial coefficients (high and low temperature)
        for all species.

        RETURNS:
        --------
        NASA_poly_coeffs: dict
            Dictionary of the form {species name : dict of NASA polynomial coefficients}
            where the inner dictionary has keys 'high' (with value, array of coeffs in high T),
            'low' (with value, array of coeffs in low T), and 'Tmid' (with value, float of middle T)
        """
        NASA_poly_coefs = {}
        for species_name in self.species:
            coef = {}
            coef['Tmid'] = self.get_Tmid(species_name)
            coef['low'] = self.get_coeffs(species_name, 'low')
            coef['high'] = self.get_coeffs(species_name, 'high')
            NASA_poly_coefs[species_name] = coef
        self.NASA_poly_coefs = NASA_poly_coefs
        return self.NASA_poly_coefs

    def get_reaction_list(self):
        """Appends a Reaction/Reaction-inherited object to each reaction
        described by input xml file to reaction_list.
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
                    rxn = Reactions.IrreversibleReaction(rxn_type, is_reversible, rxn_equation,
                                           species, rate_coeffs_components,
                                           reactant_stoich_coeffs, product_stoich_coeffs)
                    self.reaction_list.append(rxn)

                # REVERSIBLE elementary reaction case
                elif is_reversible == True and rxn_type == "Elementary":
                    rxn = Reactions.ReversibleReaction(rxn_type, is_reversible, rxn_equation,
                                         species, rate_coeffs_components,
                                         reactant_stoich_coeffs, product_stoich_coeffs)
                    self.reaction_list.append(rxn)

                # Unhandled reaction case
                else:
                    raise NotImplementedError("This type of reaction has not been implemented yet!")
