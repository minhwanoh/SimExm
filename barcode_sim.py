import numpy as np
from src.load import load_gt
from src.labeling import label
from src.optics import resolve
from src.output import save
import os
#Parameters

volume_dim = (10, 1000, 1000)
voxel_dim = (40, 8, 8)
offset = (0,0,0)
image_path = "/home/jeremy/emdata/kasthuri/ac3"
regions = {'synapse': {'region_path': "/home/jeremy/emdata/kasthuri/ac3synapse"}}

fluorophores = ['ATTO390', 'ATTO425', 'ATTO488', 'ATTO550']
synapse_fluorophore = 'ATTO700'
length_barcode = 15
amplification = 50
num_barcodes_synapse = 2
num_barcodes_cytosol = 10

synapse_params = {'layer1':
                            {'fluorophore' : synapse_fluorophore,
                            'region' : 'synapse',
                            'labeling_density' : 1.0,
                            'protein_density' : 1e-4,
                            'protein_noise' : 0.2,
                            'antibody_amp' : 5.0,
                            'single_neuron' : False
                            }
                }

expansion_params = {'factor': 4}

optics_params = {'numerical_aperture' : 1.15,
                 'focal_plane_depth' : 500,
                 'objective_back_aperture' : 1.0,
                 'exposure_time' : 0.1,
                 'objective_efficiency' : 0.8,
                 'detector_efficiency' : 0.6,
                 'objective_factor' : 40.0,
                 'pixel_size' : 6500,
                 'baseline_noise' : 0,
                 'channels': {'channel1': { 'laser_wavelength' : 390,
                                            'laser_power' : 50.0,
                                            'laser_percentage' : 0.25,
                                            'laser_filter' : [380, 400] },
                              'channel2': { 'laser_wavelength' : 425,
                                            'laser_power' : 50.0,
                                            'laser_percentage' : 0.25,
                                            'laser_filter' : [415, 435] },
                              'channel3': { 'laser_wavelength' : 488,
                                            'laser_power' : 50.0,
                                            'laser_percentage' : 0.25,
                                            'laser_filter' : [478, 498] },
                              'channel4': { 'laser_wavelength' : 550,
                                            'laser_power' : 50.0,
                                            'laser_percentage' : 0.25,
                                            'laser_filter' : [540, 560] },
                              'channel5': { 'laser_wavelength' : 700,
                                            'laser_power' : 50.0,
                                            'laser_percentage' : 0.25,
                                            'laser_filter' : [710, 730] }
                            }
                }

print "Loading data..."
gt_dataset = load_gt(image_path, offset, volume_dim, regions)

print "Labeling..."
(d, w, h) = volume_dim
volume = np.zeros((length_barcode, len(fluorophores) + 1, d, w, h), np.uint32)

#First label all synapses
synapse_volume, cells = label(gt_dataset, volume_dim, voxel_dim, synapse_params)
for i in range(length_barcode):
    volume[i,len(fluorophores), :, :, :] = synapse_volume[synapse_fluorophore]

#Next, put barocodes at the synapses
barcodes = {cell_id: np.random.randint(0, len(fluorophores), size=length_barcode) for cell_id in gt_dataset}
for cell_id in gt_dataset:
    reg = np.array(gt_dataset[cell_id]['synapse'])
    if len(reg) == 0: continue
    #barcode synapse
    num_barcodes = np.random.poisson(num_barcodes_synapse)
    indices = np.random.randint(0, len(reg), size=num_barcodes)
    voxels = reg[indices, :]
    for i in range(length_barcode):
        fluo_index = barcodes[cell_id][i]
        for (z, x, y) in voxels:
            ab_std = 50.0 / (np.array(voxel_dim) * expansion_params['factor'])#keeps it below 200nm on each side
            Z = np.clip(np.round(np.random.normal(z, ab_std[0], amplification)).astype(np.uint32), 0, d-1)
            X = np.clip(np.round(np.random.normal(x, ab_std[0], amplification)).astype(np.uint32), 0, d-1)
            Y = np.clip(np.round(np.random.normal(y, ab_std[0], amplification)).astype(np.uint32), 0, d-1)
            #Get random subset
            for j in range(len(Z)): 
                volume[i, fluo_index, Z[j], X[j], Y[j]] += 1
    #barcode cytosol
    reg = np.array(gt_dataset[cell_id]['cytosol'])
    num_barcodes = np.random.poisson(num_barcodes_cytosol)
    indices = np.random.randint(0, len(reg), size=num_barcodes)
    voxels = reg[indices, :]
    for i in range(length_barcode):
        fluo_index = barcodes[cell_id][i]
        for (z, x, y) in voxels:
            ab_std = 50.0 / (np.array(voxel_dim) * expansion_params['factor'])#keeps it below 200nm on each side
            Z = np.clip(np.round(np.random.normal(z, ab_std[0], amplification)).astype(np.uint32), 0, d-1)
            X = np.clip(np.round(np.random.normal(x, ab_std[0], amplification)).astype(np.uint32), 0, d-1)
            Y = np.clip(np.round(np.random.normal(y, ab_std[0], amplification)).astype(np.uint32), 0, d-1)
            #Get random subset
            for j in range(len(Z)): 
                volume[i, fluo_index, Z[j], X[j], Y[j]] += 1

#Resolve
path = '/home/jeremy/barcode_sim/'
if not os.path.isdir(path): os.mkdir(path)
print "Imaging..."
for i in range(length_barcode):
    labeled_volumes = {fluorophores[j]: volume[i, j, :, : , :] for j in range(len(fluorophores))}
    volumes = resolve(labeled_volumes, volume_dim, voxel_dim, expansion_params, optics_params)
    save(volumes, path, 'stack_' + str(i), 'splitted', 'tiff')

print "Done!"