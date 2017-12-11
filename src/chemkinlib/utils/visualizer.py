
from graphviz import Digraph
import random
import imageio
imageio.plugins.ffmpeg.download()
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

    def __init__(self, target, obj, integrate=False, time=None, cluster=False):
        #Get unique species of the system from other class function
        self.unique_species = [i[0] for i in obj.involved_species.items()]
        
        #Get system reaction types from other class function
        self.types = [i.is_reversible for i in obj.reaction_list]
        
        self.reactants = []
        #Get reactant species of the system from other class function
        for i in obj.reaction_list:
            temp_reaction = []
            for j in i.reactant_stoich_coeffs.items():
                temp_reaction.append(j[0])
            self.reactants.append(temp_reaction)
            
        self.products = []
        #Get product species of the system from other class function
        for i in obj.reaction_list:
            temp_reaction = []
            for j in i.product_stoich_coeffs.items():
                temp_reaction.append(j[0])
            self.products.append(temp_reaction)

        #Check if No# of Reactions consistent with No# of Specie Lists
        if len(self.reactants)!=len(self.products) or len(self.products)!=len(self.types):
            raise ValueError("No# of reaction system elements must be consistent.")
        
        #Get reactant concentrations of the system from other class function
        self.reactant_concentrations = obj.vis_concentrations
        
        #Check if Reactant Concentrations are Positive
        if sum([1 if i[1]<0 else 0 for i in self.reactant_concentrations.items()])!=0:
            raise ValueError("Specie Concentrations must be positive.")

        self.max_node_size = 5
        self.arrow_max_width = 5
        #If integrate flag set, get product concentrations and reaction rates at 'time', else constant defined by 
        #user and final reaction rates.
        if integrate==True:
            temp_conc = obj.step(time)[1]
            self.product_concentrations = dict([(i,temp_conc[ind]) for ind, i in enumerate(self.unique_species)])
            temp_rates = obj.get_reaction_rate()
            self.reaction_rates = dict([(i,temp_rates[ind]) for ind, i in enumerate(self.unique_species)])
        else:
            self.product_concentrations = dict([(i,self.max_node_size) for ind, i in enumerate(self.unique_species)])
            temp_rates = obj.get_reaction_rate()
            self.reaction_rates = dict([(i,temp_rates[ind]) for ind, i in enumerate(self.unique_species)])
        """
        #Check if Reactant Concentrations are Positive
        if sum([1 if i[1]<0 else 0 for i in self.product_concentrations.items()])!=0:
            raise ValueError("Specie Concentrations must be positive.")
        """
        self.fitted = False
        self.connected = False
        self.connections = []
        if cluster :
            self.cluster = True
            self.graph = Digraph(target, format='png')
            self.graph.attr('node', shape='doublecircle')
            self.graph.attr(label='Reaction Path Diagram')
            #self.graph.attr(size='20,20!')
            self.color_index = self.initialize_color_index()
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
            #self.prod_graph.attr(size='20,20!')
            self.color_index = self.initialize_color_index()
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
        
        self.fitted = True
        
    def connect(self, graphics_dict=None,
                                        size=5, separate=False):
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
        #Check if graph connected
        if self.fitted == False:
            raise AttributeError("Please call fit() method first.")

        if graphics_dict == None:
            raise AttributeError("Graphics dictionary not passed.")
            
        #Check if graphics dictionary is in readable form    
        if sum([0 if (i[1]==True or i[1]==False) else 1 for i in graphics_dict.items()])!=0:
            raise ValueError("Graphics Dictionary must contain only True (1) or False (0) values.")
        
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
            
        self.connected = True
            
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
        max_conc_reac = max(reac_conc.items(), key=lambda x: x[1])[1]
        max_conc_prod = max(prod_conc.items(), key=lambda x: x[1])[1]
        #Check if graph needs to be separated
        if separate:
            #Define Reactant Cluster
            with self.graph.subgraph(name='cluster_reactant') as c:
                c.attr(color=reac_color)
                c.attr(label='Reactants')
                for index, specie in enumerate(self.unique_species):
                    temp_size = str((reac_conc[specie]/max_conc_reac)*self.max_node_size)
                    if graphics_dict['node_color']==True:
                        c.node(specie+self.tag_reactant, **{'width':temp_size, 'height':temp_size}, color=reac_color)
                    else:
                        c.node(specie+self.tag_reactant, **{'width':temp_size, 'height':temp_size})
            #Define Product Cluster
            with self.graph.subgraph(name='cluster_product') as c:
                c.attr(color=prod_color)
                c.attr(label='Products')
                for index, specie in enumerate(self.unique_species):
                    temp_size = str((prod_conc[specie]/max_conc_prod)*self.max_node_size)
                    if graphics_dict['node_color']==True:
                        c.node(specie+self.tag_product, **{'width':temp_size, 'height':temp_size}, color=prod_color)
                    else:
                        c.node(specie+self.tag_product, **{'width':temp_size, 'height':temp_size})
        else:
            #Define Single Cluster
            for index, specie in enumerate(self.unique_species):
                temp_size = str((prod_conc[specie]/max_conc_prod)*self.max_node_size)
                if graphics_dict['node_color']==True:
                    self.graph.node(specie, **{'width':temp_size, 'height':temp_size}, color=reac_color)
                else:
                    self.graph.node(specie, **{'width':temp_size, 'height':temp_size})
        
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
        max_conc_reac = max(reac_conc.items(), key=lambda x: x[1])[1]
        max_conc_prod = max(prod_conc.items(), key=lambda x: x[1])[1]

        #Check if graph needs to be separated for reactants and products
        if separate:
            for index, specie in enumerate(self.unique_species):
                if graphics_dict['node_color']==True:
                    temp_size = str((reac_conc[specie]/max_conc_reac)*self.max_node_size)
                    self.reac_graph.node(specie+self.tag_reactant, **{'width':temp_size, 'height':temp_size}, color=reac_color)
                    temp_size = str((prod_conc[specie]/max_conc_prod)*self.max_node_size)
                    self.prod_graph.node(specie+self.tag_product, **{'width':temp_size, 'height':temp_size}, color=prod_color)
                else:
                    temp_size = str((reac_conc[specie]/max_conc_reac)*self.max_node_size)
                    self.reac_graph.node(specie+self.tag_reactant, **{'width':temp_size, 'height':temp_size})
                    temp_size = str((prod_conc[specie]/max_conc_prod)*self.max_node_size)
                    self.prod_graph.node(specie+self.tag_product, **{'width':temp_size, 'height':temp_size})
        else:
            for index, specie in enumerate(self.unique_species):
                temp_size = str((prod_conc[specie]/max_conc_prod)*self.max_node_size)
                if graphics_dict['node_color']==True:
                    self.reac_graph.node(specie, **{'width':temp_size, 'height':temp_size}, color=reac_color)
                else:
                    self.reac_graph.node(specie, **{'width':temp_size, 'height':temp_size})
        
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
        #Initiate final graphics dictionary
        graphics = {}
        
        #Get connection specific graphics
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
        
        #Check for Reversible        
        if connection[2]==True:
            graphics['dir']= 'both'
            
        return graphics
    
    def initialize_color_index(self):
        """
        Helper method to initialize different colors for each reaction index.
        """
        #Initialize color dictionary for edges and set random state
        color_dict = {}
        rstate = np.random.RandomState(9000)
        random_func = lambda: rstate.randint(0,255)
        
        #Get a number of colors randomly from hexadecimal color representation
        for i in range(25):
            color = '#%02X%02X%02X' % (random_func(),random_func(),random_func())
            color_dict[i] = color
        return color_dict
                                            
    def plot(self, test=False):
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
        #Check if graph connected
        if self.connected == False:
            raise AttributeError("Please call connect() method first.")
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
 
        #Check if image list is empty
        if img_list==[]:
            raise ValueError("Image list empty!")
            
        #Create and concatenate image clips
        clips = [ImageClip(img).set_duration(0.5) for img in img_list]
        concat_clip = concatenate_videoclips(clips, method="compose")
        concat_clip.write_videofile(target+".mp4", fps=24)
