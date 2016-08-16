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
# Download a sample database here: https://www.dropbox.com/sh/u03wsyxl12oh4no/AAD0QhlCSGvdl7QCwH0MkXJta?dl=0
# Note that the sample database contains slices 2000 to 3000 so don't use slices outside of that range

path = "/Users/Jeremy/Documents/janelia2000_3000.db" #Put the path to the database file you downloaded
ds = DatabaseSession(path) # Connects to the Database at this path, creates a new one if it can't find it

#Dataset id to use
selected_gt_id = 1
xy_bounds = [1000, 1500, 1000, 1500] #For Janelia : [0-2000, 0-2000]
z_slice = 2300 #For Janelia (downloaded above): [2000-3000]
number_of_slices = 10
expansion_factor = 20
gt_full = 1 # 1 for filled ground_truth, 0 for membrane ground_truth only


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
labelingUnit = BrainBowUnit(fluor_noise = 0.01, labeling_density = 0.1, protein_density = 0.1)

#Input optics parameteres
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

#Input params here
opticalUnit = Confocal(ds, baseline = [30, 30, 30], numerical_aperture = 1.1)


#RUN, don't modify this line
images, ground_truth = run_simulation(ds, selected_gt_id, xy_bounds, z_slice, number_of_slices, expansion_factor, labelingUnit, opticalUnit, gt_full)

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
