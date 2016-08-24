'''
run_sim.py

An example SimExm script. Guides through data loading + simulation
'''

#Indicate the path to SimExm
path_to_SimExm = '~/SimExm'
import sys; sys.path.insert(0, path_to_SimExm)

import os
from src.database.lib import load_to_array
from src.database.dataset import Dataset
from src.database.sim import SimStack, SimParams



#############################   DATA SETUP    ###############################



#Paths to db and original images
original_images_path = "~/connectomics_data/ground_truth_original_format/Janelia/images"
#add other region stacks here. See database.lib.load_image_stack
#synapses_images_path = ..
#axon_images_path = ..

#Parameters
voxel_dim = (8, 8, 8) #Jenalia dataset
bounds_wanted = (500, 500, 100) #width, height, depth of output stack
offset = (1000, 1000, 3000) # required ground truth data is loaded around this location

#Path to temporary db file for fast data loading (make sure to end with .h5)
db_path = "~/allo.h5"



###############################   SIM SETUP   ###############################

 #if True, the sim ground truth output is membrane only, otherwise cells are filled in the ground truth
gt_membrane_only = True

#Create the ground truth dataset object
gt_dataset = Dataset()

#Now create labeling unit, and give it a Dataset object to interact with ground truth data
labeling_unit = BrainbowUnit(gt_dataset)


#Now create the expansion unit, which takes the dataset as input
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

optical_unit = ConfocalUnit(num_channels = 3)
#Compute additional parameters, using the voxel_dim, the expansion factor and the wanted output bounds
bounds_required = optical_unit.compute_parameters(voxel_dim, expansion_unit.get_expansion_factor(), bounds_wanted)


###############################  LOAD DATA  ###############################

#Now if the dataset is created for the first time, we need to load in the data from file
if os.path.exists(db_path):
	#If file already exists, load from it
	gt_dataset.load_from_db(db_path)
else:
	#We can now start to load the ground truth data
	im_data = load_to_array(original_images_path, bounds_required)
	gt_dataset.load_from_image_stack(im_data)
	gt_dataset.save_to_db(db_path)



###############################   RUN SIM   ###############################



#We are ready to start the  simulation.
#You can now call label_cells as many times as you want with different parameters each time.
"""
Default:

region_type = 'full' #The region to label, could be synapses, dendrites, etc.. depending on what was loaded in the ground truth
labeling_density = 0.25
protein density = 0.80 # 1 means 1 per voxel on average
antibody_amplification_factor = 5
fluors = ['ATTO488', 'ATTO550', 'ATTO647N'] #can use any of  [Alexa350, Alexa790, ATTO390, ATTO425, ATTO430LS, ATTO465, ATTO488, ATTO490LS, ATTO550, ATTO647N, ATTO700]
fluor_noise = 0.00001
membrane_only = True
single_neuron = False
"""

labelingUnit.label_cells(region_type = 'full', fluors = ['ATTO488'], labeling_density = 1, membrane_only = True)
labelingUnit.label_cells(region_type = 'full', fluors = ['ATTO550'], labeling_density = 0.1, membrane_only = False)#Repeat this as many time as you'd like

fluo_volume = labelingUnit.get_labeled_volume()
fluo_gt = labelingUnit.get_ground_gruth(membrane_only = gt_membrane_only)

#Now pass to expansion unit

expanded_volume = expansion_unit.expand_volume(fluo_volume)
expanded_gt = expansion_unit.expand_ground_truth(fluo_gt)
expanded_dim = bounds_required * expansion_unit.get_expansion_factor()

#Now pass to optical unit

sim_volume = optical_unit.resolve_volume(expanded_volume, expanded_dim)
sim_gt = optical_unit.resolve_ground_truth(expanded_gt, expanded_dim)


#Create a SimStack and a SimParam object for easy storage and visualization
sim_params = SimParams(bounds_wanted, gt_dataset.get_parameters(), labeling_unit.get_parameters(), expansion_unit.get_parameters(), optical_unit.get_parameters())
sim_stack = SimStack(sim_volume, sim_gt, sim_param)

###############################   SAVE   ###############################

sim_name = "simulation1"
dest = "~/"

#Save or view the sim_stack, can save as tiff, gif, image_sequence (see src.database.models.sim.SimStack)
sim_stack.save_as_tiff(dest, sim_name)





