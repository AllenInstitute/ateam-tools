import os
import numpy as np

from bmtk.builder.networks import NetworkBuilder


# Step 1: Create a v1 mock network of 14 cells (nodes) with across 7 different cell "types"
net = NetworkBuilder("v1")
net.add_nodes(N=1,  # specifiy the number of cells belong to said group.
              pop_name='SST', location='L4', ei='i',  # pop_name, location, and ei are optional parameters that help's identifies properties of the cells. The modeler can choose whatever key-value pairs as they deem appropiate.
              positions=[(2.753, -464.868, -161.705)],  # The following properties we are passing in lists
              tuning_angle=[0.0],                 #  values to each individual cell
	      rotation_angle_xaxis=0.023094,
              rotation_angle_yaxis=0,
              rotation_angle_zaxis=-1.25763, # Note that the y-axis rotation is differnt for each cell (ie. given a list of size N), but with z-axis rotation all cells have the same value
              model_type='biophysical',  # The type of cell we are using
              model_template='ctdb:Biophys1.hoc',  # Tells the simulator that when building cells models use a hoc_template specially created for parsing Allen Cell-types file models. Value would be different if we were using NeuronML or different model files
              model_processing='aibs_allactive_ani_directed',  # further instructions for how to processes a cell model. In this case aibs_perisomatic is a built-in directive to cut the axon in a specific way
              dynamics_params='optim_param_313862167.json',  # Name of file (downloaded from Allen Cell-Types) used to set model parameters and channels
              morphology='Sst-IRES-Cre_Ai14-167638.05.02.01_501934125_m.swc'),  # Name of morphology file downloaded


# Step 2: We want to connect our network. Just like how we have node-types concept we group our connections into
# "edge-types" that share rules and properties
net.add_edges(source={'ei': 'i'}, target={'pop_name': 'SST'},
              connection_rule=5,
              syn_weight=6.4e-05,
              weight_function='wmax',
              distance_range=[0.0,1e+20],
              target_sections=['somatic', 'basal'],
              delay=2.0,
              dynamics_params='GABA_InhToInh.json',
              model_template='exp2syn')
net.build()
net.save(output_dir='network')


