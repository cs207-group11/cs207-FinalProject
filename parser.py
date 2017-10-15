
"""Module with parser for xml file"""

import numpy as np
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
        try:
            self.xml_filename = xml_filename
        except IOError as err:
            raise IOError("Reaction (xml) file not found!")

        self.reaction_list = []

    def __call__(self):
        pass


