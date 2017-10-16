import xml.etree.ElementTree as ET
import Reaction

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
        try:
            self.xml_filename = xml_filename
        except IOError as err:
            raise IOError("Reaction (xml) file not found!")

        self.reaction_list = []

    def __call__(self):
        """Parser all information for all reactions
        INPUTS
        ======
        self : contains xml file name: self.xml_filename 
                which is our target file to parser
        """
        
        tree = ET.parse(self.xml_filename)
        rxns = tree.getroot()
        # species
        phase = rxns.findall('phase')
        species_list = []
        for phase in rxns.findall('phase'):
            species = phase.find('speciesArray').text
        self.species = species.split()
        
        # rxns
        for reactionData in rxns.findall('reactionData'):
            for reaction in reactionData.findall('reaction'):
                # get id
                Id = reaction.get('id')
                # get is_reversible
                if reaction.get('reversible') == "yes":
                    is_reversible = True
                else:
                    is_reversible = False
                # type
                rxn_type = reaction.get('type')
                # rxn_equation
                rxn_equation = reaction.find('equation').text
                # reaction_coef
                for coef in reaction.findall('rateCoeff'):
                    for arr in coef.findall('Arrhenius'):
                        A = float(arr.find('A').text)
                        b = float(arr.find('b').text)
                        E = float(arr.find('E').text)
                rate_coeffs_components = {
                    "A":A,
                    "b":b,
                    "E":E
                }
                # reactant_stoich_coeffs
                reactant_stoich_coeffs = {}
                for reactant in reaction.find('reactants').text.split():
                    key = reactant.split(":")[0]
                    value = reactant.split(":")[1]
                    reactant_stoich_coeffs[key] = value
                # product_stoich_coeffs
                product_stoich_coeffs = {}
                for product in reaction.find('products').text.split():
                    key = product.split(":")[0]
                    value = product.split(":")[1]
                    product_stoich_coeffs[key] = value
                rxn = Reaction(rxn_type, is_reversible,
                 rxn_equation, rate_coeffs_components,
                 reactant_stoich_coeffs, product_stoich_coeffs)
                self.reaction_list.append(rxn)
        pass