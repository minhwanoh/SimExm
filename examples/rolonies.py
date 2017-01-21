import os
import sys
#Indicate the path to SimExm here
sys.path.append('/path/to/SimExm')
import numpy as np
from src.load import load_gt
from src.labeling import label
from src.optics import resolve
from src.output import save
import json

#####################################
#----  Ground turth parameters  ----#
#####################################

#Volume dimension (z, x, y) in # voxels to load as ground truth
volume_dim = (375, 710, 1100)
#Voxel dimension (z, x, y) in nanometers
voxel_dim = (40, 32, 32)
#Where to start loading the data, defaults to (0,0,0)
offset = (0,0,0)
#Path to main ground truth segmentation
image_path = "/path/to/main/gt/segmentation"
#Dict: {region -> {'region_path' -> path/to/synapse/annotation}}
regions = {'synapse': {'region_path': "/path/to/synapse/annotation"}}
#Output path
out_path = "/output/path"

#####################################
#------  Labeling parameters  ------#
#####################################

#4 orthogonal fluorophores for barcodes
fluorophores = ['ATTO390', 'ATTO425', 'ATTO488', 'ATTO550']
#Fluorophroes to use for cytosol and synapse stains
cytosol_fluorophore = 'ATTO647N'
synapse_fluorophore = 'ATTO700'
#Number of imaging rounds
length_barcode = 15
#Number of fluorescent protein per barcode
amplification = 50
#Barcode density at the synapses: can play around with this. 2.0 every 500 voxels seem to work well on Kasthuri
barcode_density_synapse = 2.0 / 500
#Barcode density in the cytosol: can play around with this. Tried 1.0 / 100 and 1.0 / 10.
barcode_density_cytosol = barcode_density_synapse / 100.0

#Labeling density: proportion of neurons labeled
#Protein density: 1 is 1 per nm^3
#Protein noise: don't use
cytosol_params = {'layer1':
                            {'fluorophore' : cytosol_fluorophore,
                            'region' : 'cytosol',
                            'labeling_density' : 0.2,
                            'protein_density' : 1e-3,
                            'protein_noise' : 0,
                            'antibody_amp' : 5.0,
                            'single_neuron' : False
                            }
                }
synapse_params = {'layer1':
                            {'fluorophore' : synapse_fluorophore,
                            'region' : 'synapse',
                            'labeling_density' : 1.0,
                            'protein_density' : 1e-3,
                            'protein_noise' : 0,
                            'antibody_amp' : 5.0,
                            'single_neuron' : False
                            }
                }

#####################################
#-----  Expansion parameters  ------#
#####################################

expansion_params = {'factor': 4}

#####################################
#------  Optics parameters  ------#
#####################################

#See README.md for details on the following parameters.

optics_params = {'type': 'confocal',
                 'numerical_aperture' : 1.0,
                 'refractory_index' : 1.33,
                 'focal_plane_depth' : 400,
                 'objective_back_aperture' : 1.0,
                 'exposure_time' : 0.1,
                 'objective_efficiency' : 0.8,
                 'detector_efficiency' : 0.6,
                 'objective_factor' : 40.0,
                 'pixel_size' : 6500,
                 'pinhole_radius': 50.0,
                 'baseline_noise' : 0,
                 'channels': {'channel1': { 'laser_wavelength' : 390,
                                            'laser_power' : 50.0,
                                            'laser_percentage' : 0.25,
                                            'laser_filter' : [290, 490] },
                              'channel2': { 'laser_wavelength' : 425,
                                            'laser_power' : 50.0,
                                            'laser_percentage' : 0.25,
                                            'laser_filter' : [325, 525] },
                              'channel3': { 'laser_wavelength' : 488,
                                            'laser_power' : 50.0,
                                            'laser_percentage' : 0.25,
                                            'laser_filter' : [388, 588] },
                              'channel4': { 'laser_wavelength' : 550,
                                            'laser_power' : 50.0,
                                            'laser_percentage' : 0.25,
                                            'laser_filter' : [450, 650] },
                              'channel5': { 'laser_wavelength' : 647,
                                            'laser_power' : 50.0,
                                            'laser_percentage' : 0.25,
                                            'laser_filter' : [547, 700] },
                              'channel6': { 'laser_wavelength' : 700,
                                            'laser_power' : 50.0,
                                            'laser_percentage' : 0.25,
                                            'laser_filter' : [600, 800] }
                            }
                }

###############################
##### DO NOT MODIFY BELOW #####
###############################

#Make sure path exists
if not os.path.isdir(out_path): os.mkdir(out_path)

print "Loading data..."
gt_dataset = load_gt(image_path, offset, volume_dim, 'image sequence', 'merged', False, regions)

print "Labeling cytosol..."
#First label a subset of cells with a cytosolic stain
cytosol_volume, cells = label(gt_dataset, volume_dim, voxel_dim, cytosol_params)
cytosol_volume = cytosol_volume[cytosol_fluorophore]
#First label all synapses
print "Labeling synapses..."
synapse_volume, cells = label(gt_dataset, volume_dim, voxel_dim, synapse_params)
synapse_volume = synapse_volume[synapse_fluorophore]
cells = cells[synapse_fluorophore]

#Next, put barocodes at the synapses
barcodes = {cell_id: np.random.randint(0, len(fluorophores), size=length_barcode) for cell_id in cells}
(d, w, h) = volume_dim

print "Labeling barcodes and imaging..."
#Create barcode locations
barcode_locs = {}
for cell_id in cells:
    barcode_locs[cell_id] = {}
    for region, barcode_density in [('synapse', barcode_density_synapse), ('cytosol', barcode_density_cytosol)]:
        reg = np.array(gt_dataset[cell_id][region])
        if len(reg) == 0: continue
        #barcode synapse
        barcode_mean = int(barcode_density * len(reg))
        num_barcodes = np.random.poisson(barcode_mean)
        #Select which voxels will provide a center for a rolonie
        indices = np.random.randint(0, len(reg), size=num_barcodes)
        barcode_locs[cell_id][region] = reg[indices, :]

#Now take images
for i in range(length_barcode):
    print "Imaging {}".format(i)
    volumes = {fluo: np.zeros(volume_dim, np.uint32) for fluo in fluorophores}
    #Add the synapse volume to all
    volumes[synapse_fluorophore] = synapse_volume
    volumes[cytosol_fluorophore] = cytosol_volume
    for cell_id in cells:
        for region in ['synapse', 'cytosol']:
            #Get barcode locations
            barcode_locations = barcode_locs[cell_id][region]
            fluo_index = barcodes[cell_id][i]
            volume = volumes[fluorophores[fluo_index]]
            for (z, x, y) in barcode_locations:
                #keeps it below 200nm on each side, such that each barcode is approximatly a 400nm dot in post-expansion space.
                ab_std = 100.0 / (np.array(voxel_dim) * expansion_params['factor'])
                #Dsitribute fluorophore proteins in a 200nm radius, and bound to the size of the volume
                Z = np.clip(np.round(np.random.normal(z, ab_std[0], amplification)).astype(np.uint32), 0, d-1)
                X = np.clip(np.round(np.random.normal(x, ab_std[1], amplification)).astype(np.uint32), 0, w-1)
                Y = np.clip(np.round(np.random.normal(y, ab_std[2], amplification)).astype(np.uint32), 0, h-1)
                #Populate volume
                np.add.at(volume, (Z, X, Y), 1)
    #Resolve
    out = resolve(volumes, volume_dim, voxel_dim, expansion_params, optics_params)
    save(out, out_path, 'stack_' + str(i), 'splitted', 'tiff')

# Output as dict:
# {
# cell_barcodes : { cell_id -> barcode }
# barcode_locations : {cell_id -> region -> locations}
# }
# Here region refers to synapse to cytosol

if out_path[-1] != '/': out_path += '/'
out_dict = {
            'cell_barcodes': {str(k): v.tolist() for k, v in barcodes.items()},
            'barcode_locations': {str(k): {k2: v2.tolist() for k2, v2 in v.items()} for k,v in barcode_locs.items()}
           }
with open(out_path + 'barcodes.json', 'w') as f:
    json.dump(out_dict, f)

print "Done!"