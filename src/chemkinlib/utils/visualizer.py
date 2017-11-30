
from graphviz import Digraph
import random
from moviepy.editor import *
import numpy as np

class ReactionPathDiagram():
    """
    Initializes the values of concentations, rates and species of a system of elementary reactions. 
    
    INPUTS
    =======
    target - Target file location of generated graph
    obj - Object of class Reaction/ReactionSystem which obtains the above values from prior calculations
    integrate - 

    EXAMPLE USAGE
    =========
    - ReactionPathDiagram("\results\target_file_location", obj_Reaction)
    - ReactionPathDiagram("\results\target_file_location", obj_Reaction_System)
    """
    """
    def __init__(self, target, obj, integrate=False, time=None, cluster=False):
        #Get unique species of the system from other class function
        self.unique_species = obj.get_unique_species()
        #Get system reaction types from other class function
        self.types = obj.get_reaction_types()
        #Get reactant species of the system from other class function
        self.reactants = obj.get_reaction_species()
        #Get reactant concentrations of the system from other class function
        self.reactant_concentrations = obj.get_reactant_concentrations()
        #Get product species of the system from other class function
        self.products = obj.get_product_species()
        #If integrate flag set, get product concentrations and reaction rates at 'time', else constant defined by 
        #user and final reaction rates.
        if integrate==True:
            self.product_concentrations, self.reaction_rates = self.calculate_product_concentrations(time)
        else:
            self.product_concentrations = obj.get_product_concentrations()
            self.reaction_rates = obj.get_reaction_rates()
        #Define Connections Variable
        self.connections = []
        self.reac_graph = Digraph(target, format='png')
        self.reac_graph.attr('node', shape='doublecircle')
        self.reac_graph.attr(color='lightgrey')
        self.reac_graph.attr(size='20,20!')
        self.reac_graph.attr(label='Reaction Path Diagram')
        self.prod_graph = Digraph('subgraph')
        self.prod_graph.attr(size='20,20!')
        self.color_index = self.initialize_color_index()
        self.arrow_max_width = 5
        self.tag_reactant = " | R"
        self.tag_product = " | P "
    """
        
    def __init__(self, target, unique_species, types, reac_species,
                 prod_species, reac_rates, reac_concentrations,
                 prod_concentrations, integrate, time=None, cluster=False):
        self.unique_species = unique_species
        self.types = types
        self.reactants = reac_species
        self.reactant_concentrations = reac_concentrations
        self.products = prod_species
        #If integrate flag set, get product concentrations at 'time', else constant defined by user.
        if integrate==True:
            self.product_concentrations, self.reaction_rates = self.calculate_product_concentrations(time)
        else:
            self.product_concentrations = prod_concentrations
            self.reaction_rates = reac_rates
        self.connections = []
        if cluster :
            self.cluster = True
            self.graph = Digraph(target, format='png')
            self.graph.attr('node', shape='doublecircle')
            self.graph.attr(label='Reaction Path Diagram')
            self.graph.attr(size='20,20!')
            self.color_index = self.initialize_color_index()
            self.arrow_max_width = 5
            self.tag_reactant = " | R"
            self.tag_product = " | P "
        else:
            self.cluster = False
            self.reac_graph = Digraph(target, format='png')
            self.reac_graph.attr('node', shape='doublecircle')
            self.reac_graph.attr(color='lightgrey')
            self.reac_graph.attr(size='20,20!')
            self.reac_graph.attr(label='Reaction Path Diagram')
            self.prod_graph = Digraph('subgraph')
            self.prod_graph.attr(size='20,20!')
            self.color_index = self.initialize_color_index()
            self.arrow_max_width = 5
            self.tag_reactant = " | R"
            self.tag_product = " | P "
    
    def fit(self):
        """
        Method to define graphical nodes for each unique specie at the reactant and 
        product end. For each connection, a "hex-tuple" of reactant, product, type of reaction, 
        reactant_reaction_rate and product_reaction_rate is defined.

        EXAMPLE USAGE
        =========
        *Prior*
        graph = ReactionPathDiagram(target, obj) 
        ---------
        graph.fit()
        """
        for index in range(len(self.types)):
            temp_type = self.types[index]
            temp_reactants = self.reactants[index]
            temp_products = self.products[index]
            for i in temp_reactants:
                temp_reac_rate = self.reaction_rates[i]
                for j in temp_products:
                    temp_prod_rate = self.reaction_rates[j]
                    connection = (i, j, temp_type, temp_reac_rate, temp_prod_rate, index)
                    self.connections.append(connection)
        
    def connect(self, graphics_dict={'node_color':False,'rate':True, 'arrow_size':False,
                                     'arrow_color':True,'init_con':True,'prod_con': False},
                                        size=1, separate=False):
        """
        Method to make defined connections between system node with specific graphics.
        
        INPUTS
        =======
        
        grahics_dict :
        'node_color' : If True, the nodes of each specie assume a color specific to the reactants and products
        'rate': If True, the reaction rate of the specie is displayed
        'arrow_size': If True, the thickness of the arrow is normalized for the reaction rate
        'arrow_color': If True, the colors of the arrows are different for individual reactions
        'init_con': If True, the size of the reactant nodes is set to initial concentration, else constant size
        'prod_con': If True, the size of the product nodes is set to final concentration*, else constant size 
        *integrator needs to be implemented for this feature
        
        size = 1, constant size of nodes 
        separate = If True, the reactant and product nodes for each specie are separate.
        

        EXAMPLE USAGE
        =========
        *Prior*
        graph = ReactionPathDiagram(target, obj)
        graph.fit() *Prior*
        graphics_dict = {'node_color':True,'rate':False,'arrow_size':False,
                        'arrow_color':True,'init_con':False,'prod_con': False}
        ---------------
        graph.connect(graphics_dict, time=None, size=1, separate = True)
        """
        
        #Display Product Concentration if True else constant
        if graphics_dict['prod_con']==True:
            prod_conc = self.product_concentrations
        else:
            prod_conc = dict([(i,size) for ind, i in enumerate(self.unique_species)])
        
        #Display Reactant Concentration if True else constant
        if graphics_dict['init_con']==True:
            reac_conc = self.reactant_concentrations
        else:
            reac_conc = dict([(i,size) for ind, i in enumerate(self.unique_species)])
              
        #Build Nodes
        if self.cluster:
            self.build_nodes_cluster(graphics_dict, separate, reac_conc, prod_conc, reac_color="Green", prod_color="Red")
        else:
            self.build_nodes_free(graphics_dict, separate, reac_conc, prod_conc, reac_color="Green", prod_color="Red")
        
        #Build Connections
        for connection in self.connections:
            if separate:
                org = connection[0]+self.tag_reactant
                dest = connection[1]+self.tag_product
            else:
                org = connection[0]
                dest = connection[1]
            graphics = self.get_graphics(graphics_dict, connection)
            if self.cluster:
                self.graph.edge(org, dest, **graphics)
            else:
                self.reac_graph.edge(org, dest, **graphics)
        
        #Add Product Subgraph
        if separate and not self.cluster:
            self.reac_graph.subgraph(self.prod_graph) 
            
    def build_nodes_cluster(self, graphics_dict, separate, reac_conc, prod_conc, reac_color, prod_color):
        """
        Helper method to build nodes with specific concentrations and graphics in cluster formation.
        
        INPUTS
        =======
        
        grahics_dict :
        'node_color' : If True, the nodes of each specie assume a color specific to the reactants and products
        'rate': If True, the reaction rate of the specie is displayed
        'arrow_size': If True, the thickness of the arrow is normalized for the reaction rate
        'arrow_color': If True, the colors of the arrows are different for individual reactions
        'init_con': If True, the size of the reactant nodes is set to initial concentration, else constant size
        'prod_con': If True, the size of the product nodes is set to final concentration*, else constant size 
        *integrator needs to be implemented for this feature
        
        separate = If True, the reactant and product nodes for each specie are separate.
        reac_conc = Initialized value from user.
        prod_conc = As calculated through integration.
        reac_color = "Green", pre-defined
        prod_color = "Red", pre-defined
        """
        if separate:
            
            with self.graph.subgraph(name='cluster_reactant') as c:
                c.attr(color=reac_color)
                c.attr(label='Reactants')
                for index, specie in enumerate([i[0] for i in self.unique_species]):
                    if graphics_dict['node_color']==True:
                        c.node(specie+self.tag_reactant, **{'width':str(reac_conc[specie]), 'height':str(reac_conc[specie])}, color=reac_color)
                    else:
                        c.node(specie+self.tag_reactant, **{'width':str(reac_conc[specie]), 'height':str(reac_conc[specie])})
            
            with self.graph.subgraph(name='cluster_product') as c:
                c.attr(color=prod_color)
                c.attr(label='Products')
                for index, specie in enumerate([i[0] for i in self.unique_species]):
                    if graphics_dict['node_color']==True:
                        c.node(specie+self.tag_product, **{'width':str(prod_conc[specie]), 'height':str(prod_conc[specie])}, color=prod_color)
                    else:
                        c.node(specie+self.tag_product, **{'width':str(prod_conc[specie]), 'height':str(prod_conc[specie])})
        else:
            for index, specie in enumerate([i[0] for i in self.unique_species]):
                if graphics_dict['node_color']==True:
                    self.graph.node(specie, **{'width':str(), 'height':str(reac_conc[specie])}, color=reac_color)
                else:
                    self.graph.node(specie, **{'width':str(prod_conc[specie]), 'height':str(prod_conc[specie])})
        
    def build_nodes_free(self, graphics_dict, separate, reac_conc, prod_conc, reac_color, prod_color):
        """
        Helper method to build nodes with specific concentrations and graphics, free positioning.
        
        INPUTS
        =======
        
        grahics_dict :
        'node_color' : If True, the nodes of each specie assume a color specific to the reactants and products
        'rate': If True, the reaction rate of the specie is displayed
        'arrow_size': If True, the thickness of the arrow is normalized for the reaction rate
        'arrow_color': If True, the colors of the arrows are different for individual reactions
        'init_con': If True, the size of the reactant nodes is set to initial concentration, else constant size
        'prod_con': If True, the size of the product nodes is set to final concentration*, else constant size 
        *integrator needs to be implemented for this feature
        
        separate = If True, the reactant and product nodes for each specie are separate.
        reac_conc = Initialized value from user.
        prod_conc = As calculated through integration.
        reac_color = "Green", pre-defined
        prod_color = "Red", pre-defined
        """
        if separate:
            for index, specie in enumerate([i[0] for i in self.unique_species]):
                if graphics_dict['node_color']==True:
                    self.reac_graph.node(specie+self.tag_reactant, **{'width':str(reac_conc[specie]), 'height':str(reac_conc[specie])}, color=reac_color)
                    self.prod_graph.node(specie+self.tag_product, **{'width':str(prod_conc[specie]), 'height':str(prod_conc[specie])}, color=prod_color)
                else:
                    self.reac_graph.node(specie+self.tag_reactant, **{'width':str(reac_conc[specie]), 'height':str(reac_conc[specie])})
                    self.prod_graph.node(specie+self.tag_product, **{'width':str(prod_conc[specie]), 'height':str(prod_conc[specie])})
        else:
            for index, specie in enumerate([i[0] for i in self.unique_species]):
                if graphics_dict['node_color']==True:
                    self.reac_graph.node(specie, **{'width':str(prod_conc[specie]), 'height':str(prod_conc[specie])}, color=reac_color)
                else:
                    self.reac_graph.node(specie, **{'width':str(prod_conc[specie]), 'height':str(prod_conc[specie])})
        
    def get_graphics(self, graphics_dict, connection):
        """
        Helper method to get specific graphics for each connection.
        
        INPUTS
        =======
        
        grahics_dict :
        'node_color' : If True, the nodes of each specie assume a color specific to the reactants and products
        'rate': If True, the reaction rate of the specie is displayed
        'arrow_size': If True, the thickness of the arrow is normalized for the reaction rate
        'arrow_color': If True, the colors of the arrows are different for individual reactions
        'init_con': If True, the size of the reactant nodes is set to initial concentration, else constant size
        'prod_con': If True, the size of the product nodes is set to final concentration*, else constant size 
        *integrator needs to be implemented for this feature
        
        connection = (reactant, product, reaction_type, reactant_reaction_rate, product_reaction_rate, reaction_index)
        """
        graphics = {}
        for i in graphics_dict.items():
            if i[0]=='rate' and i[1]==True:
                label = str(connection[3]) + ", " + str(connection[4])
                graphics['label'] = label
            elif i[0]=='arrow_size' and i[1]==True and connection[2]==False:
                max_rate = max(self.reaction_rates.items(), key=lambda x: x[1])[1]
                graphics['penwidth'] = str(abs(connection[3]/max_rate)*self.arrow_max_width)
            elif i[0]=='arrow_size' and i[1]==True and connection[2]==True:
                max_rate = max(self.reaction_rates.items(), key=lambda x: x[1])[1]
                graphics['penwidth'] = str(abs((connection[3]+connection[4])/(2*max_rate))*self.arrow_max_width)
            elif i[0]=='arrow_color' and i[1]==True:
                graphics['color'] = self.color_index[connection[5]]
        if connection[2]==True:
            graphics['dir']= 'both'
        return graphics
    
    def initialize_color_index(self):
        """
        Helper method to initialize different colors for each reaction index.
        """
        color_dict = {}
        rstate = np.random.RandomState(9000)
        random_func = lambda: rstate.randint(0,255)
        for i in range(25):
            color = '#%02X%02X%02X' % (random_func(),random_func(),random_func())
            color_dict[i] = color
        return color_dict
                                            
    def plot(self):
        """
        Method to display and save generated graph.
        
        EXAMPLE USAGE
        =========
        *Prior*
        graph = ReactionPathDiagram(target, obj) 
        graph.fit() 
        graphics_dict = {'node_color':True,'rate':False,'arrow_size':False,
                        'arrow_color':True,'init_con':False,'prod_con': False}
        graph.connect(graphics_dict, time=None, size=1, separate = True)
        -----------
        graph.plot()
        """
        #Display and save graph in directory
        if self.cluster:
            self.graph.view()
        else:
            self.reac_graph.view()
        
    def create_video(self, img_list, target):
        """
        Method to generate video of generated graph images.
        
        INPUTS
        ======
        img_list : List of image locations.
        target : Target location for storing generated video.
        
        EXAMPLE USAGE
        =========
        *Prior*
        graph = ReactionPathDiagram(target, obj) 
        graph.fit() 
        graphics_dict = {'node_color':True,'rate':False,'arrow_size':False,
                        'arrow_color':True,'init_con':False,'prod_con': False}
        graph.connect(graphics_dict, time=None, size=1, separate = True)
        graph.plot()
        -----------
        images = ['results/final1.gv.png','results/final2.gv.png', 'results/final3.gv.png' ]
        graph.create_video(images, "results/video")
        """
        clips = [ImageClip(img).set_duration(0.5) for img in img_list]
        concat_clip = concatenate_videoclips(clips, method="compose")
        concat_clip.write_videofile(target+".mp4", fps=24)
        
    def calculate_product_concentrations(self, time):
        #Method to calculate specie concentration at specific time interval using integrator
        #Optional Method
        pass
        
    def integrator(self, time, concentrations, reaction_rate_coeff, reaction_rates):
        #Method to implement integration
        #Optional Method
        pass