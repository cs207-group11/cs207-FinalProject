"""Module for classes and functions in ChemKinLib."""

import numbers
import numpy
import warnings
import os.path
import xml.etree.ElementTree as ET


class ReactionParser(object):
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

    def __call__(self):
        """Parses all information for all reactions in xml file
        by calling self.append_rxn()
        """
        self.append_rxn()
    
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
    
    def append_rxn(self):
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
                
                if is_reversible == False and rxn_type == "Elementary":
                    rxn = IrreversibleReaction(rxn_type, is_reversible, rxn_equation,
                                           species, rate_coeffs_components,
                                           reactant_stoich_coeffs, product_stoich_coeffs)
                    self.reaction_list.append(rxn)

                elif is_reversible == True and rxn_type == "Elementary":
                    rxn = ReversibleReaction(rxn_type, is_reversible, rxn_equation,
                                         species, rate_coeffs_components,
                                         reactant_stoich_coeffs, product_stoich_coeffs)
                    self.reaction_list.append(rxn)
                else:
                    raise NotImplementedError("This type of reaction has not been implemented yet!")


class ReactionSystem(object):
    """Class for a system of reactions extracted from an xml file"""
    def __init__(self, reaction_list):
        """Initializes ReactionSystem.
        
        INPUTS
        ------
        reaction_list : list[Reaction] or list[Reaction-like]
            list of Reaction or Reaction-inherited objects
        """
        self.reaction_list = reaction_list

    def set_temperature(self, T):
        """Sets temperature of the reaction system.

        INPUTS
        ------
        T : float
            Temperature of reaction

        NOTES
        -----
        POST:
            - Updates self.temperature
            - Raises ValueError if inputed temperature is non-positive
        """
        if T <= 0:
            raise ValueError("Temperature has to be a positive value!")
        
        for rxnObj in self.reaction_list:
            rxnObj.set_temperature(T)

    def set_concentrations(self, X):
        """Sets concentrations of the reactions.

        INPUTS
        ======
        X : dict
            dictionary with species and corresponding concentrations
        """
        for rxnObj in self.reaction_list:
            rxnObj.set_concentrations(X)

    # TODO: handle if not all species covered, handle if invalid concentrations entered (catch at this level)

    def get_reaction_rate(self):
        """Fetches reaction rate for each reaction.

        RETURNS
        -------
        list_rxn_rates : list[float]
            list of reaction rates of reactions in the system
        """
        reaction_rate_list = [rxnObj.compute_reaction_rate() for rxnObj in self.reaction_list]
        return numpy.sum(reaction_rate_list)

    # TODO: tests for these!




class Reaction(object):
    """Base class for an elementary (irreversible) reaction"""
    def __init__(self, rxn_type, is_reversible, rxn_equation, species_list,
                 rate_coeffs_components,
                 reactant_stoich_coeffs, product_stoich_coeffs):
        """
        Initializes Reaction
    
        INPUTS:
        -------
        rxn_type : str
            type of reaction (e.g. "Elementary")
        is_reversible : bool
            True if reaction is reversible
        rxn_equation : str
            string representation of reaction equation
        species_list : list
            list of chemical species from original xml file (useful for ordering)
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
        n_species = len(self.unique_species)
        return n_species

    def get_unique_species(self):
        """Helper function to return unique species involved
        in the reaction.
        
        RETURNS
        =======
        unique_species : list
            list of unique species in reaction
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
        ordered_concentrations = self.order_dictionaries(X)
        if (numpy.array(ordered_concentrations) < 0).any():
            raise ValueError("You cannot have negative concentrations!")
        
        self.concentrations = numpy.array(ordered_concentrations)

    def order_dictionaries(self, dictionary):
        """Helper function to order dictionaries (of concentrations,
        stoichiometric coefficients) based on ordering from species_list.
        This is to ensure a consistent ordering scheme.

        INPUTS
        ======
        dictionary : dict
            dictionary to order 

        RETURNS
        =======
        list_of_interest : list
            list of dictionary's keys in order of species_list
        """
        index_map = {v: i for i, v in enumerate(self.species_list)}
        sorted_tuple_list = sorted(dictionary.items(), key=lambda pair: index_map[pair[0]])
        list_of_interest = [element[1] for element in sorted_tuple_list]
        return list_of_interest


    def compute_reaction_rate(self, T=None):
        """Computes reaction rates of reaction.

        RETURNS
        =======
        rxn_rate_array: numpy.ndarray
            Array of reaction rates of reaction
        """

        reactant_stoich_coeffs = numpy.array(self.order_dictionaries(self.reactant_stoich_coeffs))
        product_stoich_coeffs = numpy.array(self.order_dictionaries(self.product_stoich_coeffs))
        concen_array = self.concentrations

        if (reactant_stoich_coeffs < 0).any():
            raise ValueError("Reactant stoichiometric coefficients must be positive!")
        
        if (product_stoich_coeffs < 0).any():
            raise ValueError("Product stoichiometric coefficients must be positive!")

        progress_rate = self.compute_progress_rate(T)
        nu_i = product_stoich_coeffs - reactant_stoich_coeffs

        reaction_rate_1_eq = progress_rate*nu_i
        return reaction_rate_1_eq


class IrreversibleReaction(Reaction):
    """Class for irreversible reaction"""
    def __init__(self, rxn_type, is_reversible, rxn_equation, species_list, rate_coeffs_components,
                 reactant_stoich_coeffs, product_stoich_coeffs):
        super(IrreversibleReaction, self).__init__(rxn_type, is_reversible, rxn_equation, species_list, rate_coeffs_components,
                         reactant_stoich_coeffs, product_stoich_coeffs)


    def compute_reaction_rate_coeff(self, T=None):
        """Computes reaction rate coefficients of reaction.

        RETURNS
        =======
        k : numeric type (or list of numeric type)
            Reaction rate coefficient
        """
        k = ReactionCoeff(self.rate_coeffs_components, T=self.temperature).k
        self.rxn_rate_coeff = k
        return k


    def compute_progress_rate(self, T=None):
        """Computes progress rates of reaction.

        RETURNS
        =======
        omega_array : numpy.ndarray
            Array of progress rates of reaction
        """
        reactant_stoich_coeffs = numpy.array(self.order_dictionaries(self.reactant_stoich_coeffs))
        concen_array = self.concentrations

        if len(concen_array) == 0:
            raise ValueError("You must set the concentrations first!")

        k = self.compute_reaction_rate_coeff(T)

        concen_powered_j = concen_array ** reactant_stoich_coeffs

        #TODO: move to coeffs computation
        if k < 0:
            raise ValueError("Reaction rate constants must be positive!")

        progress_rate = k * numpy.prod(concen_powered_j)


        return progress_rate


class ReversibleReaction(Reaction):
    """Class for reversible reaction"""
    def __init__(self, rxn_type, is_reversible, rxn_equation, species_list, rate_coeffs_components,
                 reactant_stoich_coeffs, product_stoich_coeffs):
        super().__init__(rxn_type, is_reversible, rxn_equation, species_list, rate_coeffs_components,
                 reactant_stoich_coeffs, product_stoich_coeffs)


    def compute_reaction_rate_coeff(self, T=None):
        """Computes reaction rate coefficients of reaction.

        RETURNS
        =======
        k : numeric type (or list of numeric type)
            Reaction rate coefficient
        """
        coeffs = ReactionCoeff(self.rate_coeffs_components, T=self.temperature)
        self.forward_rxn_rate_coeff = coeffs.k
        self.backward_rxn_rate_coeff = coeffs.k_back

        #just to conform with the crude interface
        return (self.forward_rxn_rate_coeff, self.backward_rxn_rate_coeff)



    def compute_progress_rate(self, T=None):
        """Computes progress rates of reaction.

        RETURNS
        =======
        omega_array : numpy.ndarray
            Array of progress rates of reaction
        """
        reactant_stoich_coeffs = numpy.array(self.order_dictionaries(self.reactant_stoich_coeffs))
        product_stoich_coeffs = numpy.array(self.order_dictionaries(self.product_stoich_coeffs))
        concen_array = self.concentrations

        if len(concen_array) == 0:
            raise ValueError("You must set the concentrations first!")

        #compute the forward and backward reaction coeffs
        self.compute_reaction_rate_coeff(T)
        print("reactant_stoich_coeffs", reactant_stoich_coeffs, "n_rxns", reactant_stoich_coeffs.shape)

        concen_powered_j_forward = concen_array ** reactant_stoich_coeffs
        concen_powered_j_backward = concen_array ** product_stoich_coeffs

        #compute the forward part
        progress_rate_forward = self.forward_rxn_rate_coeff * numpy.prod(concen_powered_j_forward)

        #compute the backward part
        progress_rate_backward = self.backward_rxn_rate_coeff * numpy.prod(concen_powered_j_backward)

        #code for debugging, save for future use
        #print("progress_rate_forward", progress_rate_forward, "progress_rate_backward", progress_rate_backward)
        #print("concen_powered_j_forward", concen_powered_j_forward, "concen_powered_j_backward", concen_powered_j_backward)
        #print("concen_array", concen_array, "reactant_stoich_coeffs", reactant_stoich_coeffs, "product_stoich_coeffs", product_stoich_coeffs)

        return progress_rate_forward-progress_rate_backward



class ReactionCoeff(object):
    """Class for reaction rate coefficients, or values k."""
    def __init__(self, k_parameters, T=None):
        """Initializes reaction rate coefficients.

        INPUTS
        ======
        T : int or float
            temperature of the reaction (in Kelvin)
        k_parameters : dictionary
            dictionary of parameters to compute k
        """
        self.k_parameters = k_parameters
        self.T = T
        self.k = self.get_coeff(self.k_parameters, self.T)
        #TODO: change when the difinition of get_back_coeff() changes
        self.k_back = self.get_back_coeff()

    def get_coeff(self, k_parameters, T):
        """Computes reaction rate coefficients depending on passed parameters.

        INPUTS
        ======
        T : int or float
            temperature of the reaction (in Kelvin)
        k_parameters : dictionary
            dictionary of parameters to compute k
        
        RETURNS
        =======
        k : int or float
            reaction rate coefficient of the reaction

        NOTES
        =====
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


    def get_back_coeff(self):
        #TODO: placeholder for the back_coeff
        return 1e-27


    def const(self, k):
        """Returns constant reaction rate coefficients k.

        INPUTS
        =======
        k : numeric type 
            constant reaction rate coefficient

        RETURNS
        ========
        k : numeric type
            constant reaction rate coefficients.

        NOTES
        =====
        POST:
            - Raises ValueError if k is non-positive!
        """
        if k <= 0:
            raise ValueError("Reaction rate must be positive.")
        return k

    def arr(self, A, E, T, R=8.314):
        """Returns Arrhenius reaction rate coefficients k.

        INPUTS
        ======
        A: float, strictly positive, no default value
           The Arrhenius prefactor
        E: float, no default value
           The Arrhenius parameter
        T: float, strictly positive, no default value
           Temperature T, asuuming a Kelvin scale
        R: float, default value is 8.314, cannot be changed except to convert units
           The ideal gas constant

        RETURNS
        =======
        k: Arrhenius reaction rate coefficients k,
           floats
           unless A or T is not postive
           in which case a ValueError exception is raised

        NOTES
        =====
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

        NOTES
        =====
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


# if __name__ == "__main__":

#     xml_filename = "rxns.xml"
#     parser = ReactionParser(xml_filename)
#     parser()
#     rxn1 = parser.reaction_list[0]

#     print rxn1
#     # print rxn1.species_list
#     # print rxn1.reactant_stoich_coeffs
#     # print rxn1.product_stoich_coeffs
#     # print rxn1.rate_coeffs_components

#     rxn1.set_concentrations({'H':1, 'O2':2, 'OH':0, 'O':0, 'H2O':0, 'H2':0})
#     rxn1.set_temperature(100)
#     rxn1.compute_reaction_rate_coeff()
#     omega = rxn1.compute_progress_rate()
#     rxnrate = rxn1.compute_reaction_rate()

#     print("Reaction rate coefficient (k) : {}".format(rxn1.rxn_rate_coeff))
#     print("Reaction progress rate: {}".format(omega))
#     print("Reaction rate: {}".format(rxnrate))

