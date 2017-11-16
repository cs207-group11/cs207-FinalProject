
import numbers
import numpy
import warnings
import os.path

class ReactionSystem():
    """Class for a system of reactions"""
    def __init__(self, reaction_list):
        self.reaction_list = reaction_list

    def set_temperature(self, T):
        """Sets temperature of the reaction system
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
        for r in self.reaction_list:
            r.set_temperature(T)

    def set_concentrations(self, X):
        """Sets concentrations of the reaction
        INPUTS
        ======
        X : dict
            dictionary with species and corresponding concentrations
        """
        for r in self.reaction_list:
            r.set_concentrations(X)

    def get_reaction_rate(self):
        reaction_rate_list = [r.compute_reaction_rate() for r in self.reaction_list]
        return sum(reaction_rate_list)




class Reaction():
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
        super().__init__(rxn_type, is_reversible, rxn_equation, species_list, rate_coeffs_components,
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
                 reactant_stoich_coeffs, product_stoich_coeffs,NASA_poly_coefs):
        super().__init__(rxn_type, is_reversible, rxn_equation, species_list, rate_coeffs_components,
                 reactant_stoich_coeffs, product_stoich_coeffs)
        self.NASA_poly_coefs = NASA_poly_coefs

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



class ReactionCoeff():
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

    ############################## format NASA ##########################
    def get_NASA_poly_coef(self):
        NASA = []
        for nasa in self.NASA_poly_coef:
            if self.temperature <= nasa["Tmid"]: #get the low temperature
                NASA.append(nasa["low"])
            else:
                NASA.append(nasa["high"])
        self.nasa7_coeffs = np.array(NASA)
        
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