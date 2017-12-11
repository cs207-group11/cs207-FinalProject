

"""Test module for Visualizer."""


import pytest
import os
import imageio
imageio.plugins.ffmpeg.download()
from chemkinlib.utils import Parser
from moviepy.editor import *
from chemkinlib.utils import visualizer
from chemkinlib.reactions import ReactionSystems
from chemkinlib.config import DATA_DIRECTORY
import numpy

# ======================= TESTS FOR REACTION-PATH-DIAGRAM OBJECT ====================== #

    
def test_fitted():
    """Test if No# of Reactions consistent with No# of Specie Lists"""
    xml_filename = os.path.join(DATA_DIRECTORY, "rxnset_long.xml")
    parser = Parser.ReactionParser(xml_filename)
    concentration = ({'H':1, 'H2':1, 'H2O':0, 'H2O2':1, 'HO2':1, 'O':1, "O2":1, "OH":1})
    temperature = 1000
    rxnsys = ReactionSystems.ReactionSystem(parser.reaction_list,
                        parser.NASA_poly_coefs,
                        temperature,
                        concentration)
    target = "final/results"
    graphics_dict = {'node_color':False,'rate':False,'arrow_size':False,'arrow_color':True,'init_con':True,'prod_con': True}
    obj = visualizer.ReactionPathDiagram(target, rxnsys, integrate=False, time=None, cluster=False)
    with pytest.raises(AttributeError):
        obj.connect()
        
    
def test_connected():
    """Test if No# of Reactions consistent with No# of Specie Lists"""
    xml_filename = os.path.join(DATA_DIRECTORY, "rxnset_long.xml")
    parser = Parser.ReactionParser(xml_filename)
    concentration = ({'H':1, 'H2':1, 'H2O':0, 'H2O2':1, 'HO2':1, 'O':1, "O2":1, "OH":1})
    temperature = 1000
    rxnsys = ReactionSystems.ReactionSystem(parser.reaction_list,
                        parser.NASA_poly_coefs,
                        temperature,
                        concentration)
    target = "final/results"
    graphics_dict = {'node_color':False,'rate':False,'arrow_size':False,'arrow_color':True,'init_con':True,'prod_con': True}
    obj = visualizer.ReactionPathDiagram(target, rxnsys, integrate=False, time=None, cluster=False)
    obj.fit()
    with pytest.raises(AttributeError):
        obj.plot()

def test_graphics_dict():
    """Test if No# of Reactions consistent with No# of Specie Lists"""
    xml_filename = os.path.join(DATA_DIRECTORY, "rxnset_long.xml")
    parser = Parser.ReactionParser(xml_filename)
    concentration = ({'H':1, 'H2':1, 'H2O':0, 'H2O2':1, 'HO2':1, 'O':1, "O2":1, "OH":1})
    temperature = 1000
    rxnsys = ReactionSystems.ReactionSystem(parser.reaction_list,
                        parser.NASA_poly_coefs,
                        temperature,
                        concentration)
    target = "final/results"
    graphics_dict = {'node_color':"hello",'rate':False,'arrow_size':False,'arrow_color':True,'init_con':True,'prod_con': True}
    obj = visualizer.ReactionPathDiagram(target, rxnsys, integrate=False, time=None, cluster=False)
    obj.fit()
    with pytest.raises(ValueError):
        obj.connect(graphics_dict, size=5, separate = False)

def test_graphics():
    """Test if No# of Reactions consistent with No# of Specie Lists"""
    xml_filename = os.path.join(DATA_DIRECTORY, "rxnset_long.xml")
    parser = Parser.ReactionParser(xml_filename)
    concentration = ({'H':1, 'H2':1, 'H2O':0, 'H2O2':1, 'HO2':1, 'O':1, "O2":1, "OH":1})
    temperature = 1000
    rxnsys = ReactionSystems.ReactionSystem(parser.reaction_list,
                        parser.NASA_poly_coefs,
                        temperature,
                        concentration)
    target = "final/results"
    graphics_dict = {'node_color':True,'rate':False,'arrow_size':False,'arrow_color':True,'init_con':True,'prod_con': True}
    obj = visualizer.ReactionPathDiagram(target, rxnsys, integrate=False, time=None, cluster=False)
    obj.fit()
    with pytest.raises(AttributeError):
        obj.connect(size=5, separate = False)
