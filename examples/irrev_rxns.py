
"""Example of irreversible reaction."""

import os
from chemkinlib.utils import Parser
from chemkinlib.utils import visualizer
from chemkinlib.reactions import ReactionSystems
from chemkinlib.config import DATA_DIRECTORY
import numpy

# USER INPUT: reaction (xml) file
xml_filename = os.path.join(DATA_DIRECTORY, "rxnset_long.xml")

parser = Parser.ReactionParser(xml_filename)

# USER INPUTS (concentrations and temperatures)
concentration = ({'H':1, 'H2':1, 'H2O':0, 'H2O2':1, 'HO2':1, 'O':1, "O2":1, "OH":1})
temperature = 1000


# Set up reaction system
rxnsys = ReactionSystems.ReactionSystem(parser.reaction_list,
                        parser.NASA_poly_coefs,
                        temperature,
                        concentration)

target = "outputs/"
graph_dict = {'node_color':True,'rate':False,'arrow_size':True,'arrow_color':True,'init_con':True,'prod_con': True}
#compute the concentration change with timestep

for i in range(3):
    graph = visualizer.ReactionPathDiagram(target+str(i), rxnsys, integrate=True, time=1e-15, cluster=True)
    graph.fit()
    graph.connect(graph_dict, size=5, separate = True)
    graph.plot()

imgs = [target+str(i)+".gv.png" for i in range(3)]
graph.create_video(imgs, target+'output_video')

# Compute and sort reaction rates
rxnrates_dict = rxnsys.sort_reaction_rates()

# display reaction rates by species
for k, v in rxnrates_dict.items():
    print("d[{0}]/dt : \t {1:e}".format(k, v))




