'''
run_sim.py

An example SimExm script. Guides through data loading + simulation
'''

#Indicate the path to SimExm
path_to_SimExm = '/home/jeremy/SimExm'
import sys; sys.path.insert(0, path_to_SimExm)

import os
from src.database.lib import load_to_array
from src.database.models.dataset import Dataset
from src.database.models.sim import SimStack, SimParams
from src.simulation.models.expansion import ExpansionUnit
from src.simulation.models.labeling import BrainbowUnit
from src.simulation.models.optics import ConfocalUnit
from PIL import Image
import numpy as np



#############################   DATA SETUP    ###############################



#Paths to db and original images
original_images_path = "/home/jeremy/janelia/ground_truth"
#add other region stacks here. See database.lib.load_image_stack
#synapses_images_path = ..
#axon_images_path = ..

#Parameters
voxel_dim = (8, 8, 8) #Jenalia dataset
bounds_wanted = (5, 500, 500) #depth, width, height of output stack
offset = (0, 1000, 1000) # (z, x, y) offset from which to load the images (crops the stack)

#Path to temporary db file for fast data loading (make sure to end with .hdf5)
db_path = "/home/jeremy/allo.hdf5"



###############################   SIM SETUP   ###############################

 #if True, the sim ground truth output is membrane only, otherwise cells are filled in the ground truth
gt_membrane_only = True

#Create the ground truth dataset object
gt_dataset = Dataset()

#Now create labeling unit, and give it a Dataset object to interact with ground truth data
labeling_unit = BrainbowUnit(gt_dataset)


#Now create the expansion unit
"""
Default:

expansion_factor = 20
"""
expansion_unit = ExpansionUnit(expansion_factor = 20)


#Now create optical unit and adjust parameters
#Note that setting the number of channels to a lower value automatically truncates laser_wavelenths, laser_power, etc.. at running time


"""
Default:

num_channels =  3

laser_wavelengths = [488, 565, 660]
laser_power  = [50, 50, 50]
laser_percentage   = [0.25, 0.25, 0.25]
objective_back_aperture  = 1.0
baseline = [50, 50, 50]
filters = [[450, 550], [550, 625], [625, 950]]
exposure_time  = 0.1
numerical_aperture  = 1.15
objective_efficiency  = 0.8
detector_efficiency  = 0.6
focal_plane_depth  = 500
objective_factor  = 40.0
pixel_size = 6500
"""

optical_unit = ConfocalUnit(num_channels = 3, baseline = [0,0,0], numerical_aperture = 1.1)

#Compute additional parameters, using the voxel_dim, the expansion factor and the wanted output bounds
bounds_required = optical_unit.compute_parameters(voxel_dim, expansion_unit.get_expansion_factor(), bounds_wanted)

###############################  LOAD DATA  ###############################

#Now if the dataset is created for the first time, we need to load in the data from file
if os.path.exists(db_path):
	#If file already exists, load from it
	gt_dataset.load_from_db(db_path)
else:
	#If not, we load from data path
	print "Loading data..."
	im_data = load_to_array(original_images_path, offset, bounds_required)
	gt_dataset.load_from_image_stack(voxel_dim, im_data)
	gt_dataset.save_to_db(db_path) #Save for quick reuse



###############################   RUN SIM   ###############################



#We are ready to start the  simulation.
#You can now call label_cells as many times as you want with different parameters each time.
"""
Default:

region_type = 'full' #The region to label, could be synapses, dendrites, etc.. depending on what was loaded in the ground truth
labeling_density = 0.25
protein density = 0.80 # 1 means 1 fluorophore per nm3
antibody_amplification_factor = 5
fluors = ['ATTO488', 'ATTO550', 'ATTO647N'] #can use any of  [Alexa350, Alexa790, ATTO390, ATTO425, ATTO430LS, ATTO465, ATTO488, ATTO490LS, ATTO550, ATTO647N, ATTO700]
fluor_noise = 0.00001
membrane_only = True
single_neuron = False
"""

print "Performing labeling simulation..."

#Repeat this as many time as you'd like
labeling_unit.label_cells(region_type = 'full', fluors = ['ATTO488'], labeling_density = 1, protein_density = 0.6, membrane_only = True)
labeling_unit.label_cells(region_type = 'full', fluors = ['ATTO550', 'ATTO647N'], labeling_density = 0.1, membrane_only = False)


fluo_volume = labeling_unit.get_labeled_volume()
fluo_gt = labeling_unit.get_ground_truth(membrane_only = gt_membrane_only)

#Now pass to expansion unit

print "Performing expansion simulation..."

expanded_volume = expansion_unit.expand_volume(fluo_volume)
expanded_gt = expansion_unit.expand_ground_truth(fluo_gt)
expanded_dim = [bounds_required[i] * expansion_unit.get_expansion_factor() for i in range(len(bounds_required))]

print "Performing optics simulation..."

#Now pass to optical unit
sim_volume = optical_unit.resolve_volume(expanded_volume, expanded_dim, labeling_unit.get_fluors_used())
sim_gt = optical_unit.resolve_ground_truth(expanded_gt, expanded_dim)


#Create a SimStack and a SimParam object for easy storage and visualization
sim_params = SimParams(bounds_wanted, gt_dataset.get_parameters(), labeling_unit.get_parameters(),\
	expansion_unit.get_parameters(), optical_unit.get_parameters())
sim_stack = SimStack(sim_volume, sim_gt, sim_params)

###############################   SAVE   ###############################

print "Saving results..."

sim_name = "simulation1"
dest = "/home/jeremy/"

#Save or view the sim_stack, can save as tiff, gif, image_sequence (see src.database.models.sim.SimStack)
sim_stack.save_as_tiff(dest, sim_name)

print "Done!"




