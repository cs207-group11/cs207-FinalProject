

"""Test module for Visualizer."""

from graphviz import Digraph
import random
import numpy as np
import os.path
import pytest
from chemkinlib.utils import visualizer

# ======================= TESTS FOR REACTION-PATH-DIAGRAM OBJECT ====================== #

def test_reaction_count():
    """Test if No# of Reactions consistent with No# of Specie Lists"""
    unique = ['a','b','c','d']
    types = [True, False]
    a = [['a','b'],['c'],['d']]
    b = [['c'], ['b','d']]
    conc_r = {'a':1,'b':2,'c':3,'d':0}
    conc_p = {'a':1.6,'b':3.1,'c':0.5,'d':5}
    rates = {'a':-2.5,'b':2, 'c':8.5, 'd':10.5}
    target = 'results/test'
    with pytest.raises(ValueError):
        obj = visualizer.ReactionPathDiagram(target, unique, types, a, b, rates, conc_r, conc_p, integrate=False, time=None, cluster=False)

def test_reactant_concetration():
    """Test if Reactant Concentrations are Positive"""
    unique = ['a','b','c','d']
    types = [True, False]
    a = [['a','b'],['c']]
    b = [['c'], ['b','d']]
    conc_r = {'a':-1,'b':2,'c':3,'d':0}
    conc_p = {'a':1.6,'b':3.1,'c':0.5,'d':5}
    rates = {'a':-2.5,'b':2, 'c':8.5, 'd':10.5}
    target = 'results/test'
    with pytest.raises(ValueError):
        obj = visualizer.ReactionPathDiagram(target, unique, types, a, b, rates, conc_r, conc_p, integrate=False, time=None, cluster=False)
    
def test_prod_concentration():
    """Test if Product Concentrations are Positive"""
    unique = ['a','b','c','d']
    types = [True, False]
    a = [['a','b'],['c']]
    b = [['c'], ['b','d']]
    conc_r = {'a':1,'b':2,'c':3,'d':0}
    conc_p = {'a':1.6,'b':-3.1,'c':0.5,'d':5}
    rates = {'a':-2.5,'b':2, 'c':8.5, 'd':10.5}
    target = 'results/test'
    with pytest.raises(ValueError):
        obj = visualizer.ReactionPathDiagram(target, unique, types, a, b, rates, conc_r, conc_p, integrate=False, time=None, cluster=False)
    
def test_fitted():
    """Test if graph is fitted."""
    unique = ['a','b','c','d']
    types = [True, False]
    a = [['a','b'],['c']]
    b = [['c'], ['b','d']]
    conc_r = {'a':1,'b':2,'c':3,'d':0}
    conc_p = {'a':1.6,'b':3.1,'c':0.5,'d':5}
    rates = {'a':-2.5,'b':2, 'c':8.5, 'd':10.5}
    target = 'results/test'
    obj = visualizer.ReactionPathDiagram(target, unique, types, a, b, rates, conc_r, conc_p, integrate=False, time=None, cluster=False)
    with pytest.raises(AttributeError):
        obj.connect()
        
    
def test_connected():
    """Test if graph is connected."""
    unique = ['a','b','c','d']
    types = [True, False]
    a = [['a','b'],['c']]
    b = [['c'], ['b','d']]
    conc_r = {'a':1,'b':2,'c':3,'d':0}
    conc_p = {'a':1.6,'b':3.1,'c':0.5,'d':5}
    rates = {'a':-2.5,'b':2, 'c':8.5, 'd':10.5}
    target = 'results/test'
    obj = visualizer.ReactionPathDiagram(target, unique, types, a, b, rates, conc_r, conc_p, integrate=False, time=None, cluster=False)
    obj.fit()
    with pytest.raises(AttributeError):
        obj.plot()

def test_graphics_dict():
    """Test if graphics_dict is readable."""
    unique = ['a','b','c','d']
    types = [True, False]
    a = [['a','b'],['c']]
    b = [['c'], ['b','d']]
    conc_r = {'a':1,'b':2,'c':3,'d':0}
    conc_p = {'a':1.6,'b':3.1,'c':0.5,'d':5}
    rates = {'a':-2.5,'b':2, 'c':8.5, 'd':10.5}
    target = 'results/test'
    graphics_dict = {'node_color':77,'rate':False,'arrow_size':True,'arrow_color':True,'init_con':True,'prod_con': False}
    obj = visualizer.ReactionPathDiagram(target, unique, types, a, b, rates, conc_r, conc_p, integrate=False, time=None, cluster=False)
    obj.fit()
    with pytest.raises(ValueError):
        obj.connect(graphics_dict, size=1, separate = False)
    
def test_video_image_list():
    """Test if image list for video creation is empty."""
    unique = ['a','b','c','d']
    types = [True, False]
    a = [['a','b'],['c']]
    b = [['c'], ['b','d']]
    conc_r = {'a':1,'b':2,'c':3,'d':0}
    conc_p = {'a':1.6,'b':3.1,'c':0.5,'d':5}
    rates = {'a':-2.5,'b':2, 'c':8.5, 'd':10.5}
    target = 'results/test'
    graphics_dict = {'node_color':True,'rate':False,'arrow_size':True,'arrow_color':True,'init_con':True,'prod_con': False}
    obj = visualizer.ReactionPathDiagram(target, unique, types, a, b, rates, conc_r, conc_p, integrate=False, time=None, cluster=False)
    obj.fit()
    obj.connect(graphics_dict, size=1, separate = False)
    obj.plot()
    imgs = []
    with pytest.raises(ValueError):
        obj.create_video(imgs, "results/test")
