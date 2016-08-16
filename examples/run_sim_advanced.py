#Select SimExm path
path_to_SimExm = ''

import sys; sys.path.insert(0, path_to_SimExm)

#Other imports
from database.models import DatabaseSession
from simulation.methods import run_simulation
from scipy.misc import toimage
import numpy as np
from math import ceil, sqrt
from simulation.models.labeling import BrainBowUnit, BrainbowPlusCytosolic
from simulation.models.optics import Confocal
from visualization.methods import make_gif, save_image_sequence, make_tiff_stack

#Path to database on local machine

path = "/Users/Jeremy/Documents/janelia2000_3000.db" #Put the path to the database file you downloaded
ds = DatabaseSession(path) # Connects to the Database at this path, creates a new one if it can't find it

#Dataset id to use
selected_gt_id = 1
xy_bounds = [1000, 1500, 1000, 1500]
center_z_slice = 2300
number_of_slices = 10
expansion_factor = 20

#Find number of z-slices required for simulation given number of slices wanted at output

voxel_dims = np.array(d.get_voxel_dimensions_nm(ds, gt_dataset_id), np.int)
slices_required = optical_unit.calculate_parameters(voxel_dims, number_of_slices)
slices_required = int(ceil(float(slices_required) / float(expansion_factor)))

#Computed bounds
ground_truth_bounds = [bounds[0], bounds[1], bounds[2], bounds[3], center_z_slice - slices_required / 2, center_z_slice + slices_required / 2]
volume_dim = [real_bounds[1] - real_bounds[0], real_bounds[3] - real_bounds[2], real_bounds[5] - real_bounds[4]]

#Load data
membranes = c.get_voxels_in_volume_by_cell(ds, gt_dataset_id, ground_truth_bounds, include_body = False,  include_membrane = True)
body = c.get_voxels_in_volume_by_cell(ds, gt_dataset_id, ground_truth_bounds, include_body = True,  include_membrane = False)

#Input labeling parameteres
"""
Default:

labeling_density = 0.25
protein density = 0.80 # 1 means 1 per voxel on average
antibody_amplification_factor = 5
fluors = ['ATTO488', 'ATTO550', 'ATTO647N'] #can use any of  [Alexa350, Alexa790, ATTO390, ATTO425, ATTO430LS, ATTO465, ATTO488, ATTO490LS, ATTO550, ATTO647N, ATTO700]
fluor_noise = 0.00001
"""

#Input params here
#labelingUnit = BrainBowUnit(ds, fluor_noise = 0.01, labeling_density = 0.1, protein_density = 0.1)

#2 color brainbow
labelingUnit_1 = BrainbowUnit(protein_density = 0.8, labeling_density = 0.9, fluors = ['ATTO425', 'ATTO550'],\
       antibody_amplification_factor = 5, fluor_noise = 0.001)

#The labeling method takes a dictionary with the voxels to annotate as first argument and a dictionary for the ground truth as second argument
(fluorophores_1, fluorophores_gt_1) = labelingUnit_1.perform_labeling(membranes, body, volume_dim, voxel_dims)

#Let say we want to give a single neuron a cytolosolic stain 
    
labelingUnit_2 = BrainbowUnit(protein_density = 0.05, labeling_density = 1.0, fluors = ['ATTO647N'],\
       antibody_amplification_factor = 5, fluor_noise = 0.00001)
       
all_neurons = np.unique(fluorophores_gt_1) #Select one of the neurons stained by the first labeling unit
selected_neuron = all_neurons[randint(0, all_neurons.size - 1)]
single_neuron_dict = {selected_neuron : body[selected_neuron]} #Select all of that neuron's voxels

(fluorophores_2, fluorophores_gt_2) = labelingUnit_2.perform_labeling(single_neuron_dict, single_neuron_dict, volume_dim, voxel_dims)

fluorophores = np.concatenate(fluorohores_1, fluorophores_2, 0) #Stack to 3-channel volume (2 + 1 here)
fluorophores_gt = np.add(fluorophores_gt_1, fluorophores_gt_2) # Merge the ground truth volumes 

#Expansion
print "Perform expansion..."
volume = perform_expansion(fluorophores, expansion_factor)
volume_gt = perform_expansion_gt(fluorophores_gt, expansion_factor)

#Input optics parameteres
"""
Default:

num_channels = 3
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
precision  = 1.0
"""

#Input params here
opticalUnit = Confocal(ds, baseline = [30, 30, 30],  numerical_aperture = 1.1)

print "Perform optics..."
channels = []
new_dims = [expansion_factor * i for i in dims]
for channel in range(opticalUnit.num_channels):
    fluo_volume = optical_unit.resolve_volume(volume, new_dims, labeling_unit.fluor_types_used, channel)
    channels.append(fluo_volume)

fluo_volume_gt = optical_unit.resolve_ground_truth(volume_gt, new_dims)

#fake channel for rgb mode. If you're only using two channels in the simulation.
#channels.append(np.empty_like(fluo_volume)) 

#Stack channels
images = []
ground_truth = []
for i in range(0, fluo_volume.shape[0]):
    images.append(np.stack([np.squeeze(fluorophore_volume[i,:,:]) for fluorophore_volume in channels], axis=-1))
    ground_truth.append(fluo_volume_gt[i, :, :])

#Destination path
dest = "/Users/Jeremy/Desktop/"

#To save as GIF:
#make_gif(images, path=dest, name="simulation1")
#make_gif(ground_truth, path=dest, name="simulation1_gt")

#To save as TIFF stack:
make_tiff_stack(images, path=dest, name="simulation_train", rgb=True, sixteen_bit_mode = False)
make_tiff_stack(ground_truth, path=dest, name="simulation_train_gt", rgb=False, sixteen_bit_mode = False)

#To save as image sequence:
#dest_gt = "/Users/Jeremy/Desktop/simulation3_gt"
#save_image_sequence(images, dest, rgb=True)
#save_image_sequence(images, dest_gt, rgb=False)
