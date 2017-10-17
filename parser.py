import xml.etree.ElementTree as ET
import os.path
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
        """
        if os.path.isfile(xml_filename):
            self.xml_filename = xml_filename
        else:
            raise IOError("Reaction (xml) file not found!")
                      
#         NOTES:
#         ------
#         POST:
#             - Raises IOError if inputed xml file not found
#                 except IOError as err:
#             raise IOError("Reaction (xml) file not found!")
        
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
                    # Arrhenius
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
                    # modified Arrhenius
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
                    # constant
                    if coef.find('Constant') is not None:
                        for arr in coef.findall('Constant'):
                            if arr.find('k') is None or arr.find('k').text is None:
                                raise ValueError("Didn't provide coefficient k for modified Arrhenius ")
                            else:
                                k = float(arr.find('k').text)
                        rate_coeffs_components = {
                            "k":k
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