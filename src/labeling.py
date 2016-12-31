# Provided under the MIT License (MIT)
# Copyright (c) 2016 Jeremy Wohlwend

# Permission is hereby granted, free of charge, to any person obtaining 
# a copy of this software and associated documentation files (the "Software"),
# to deal in the Software without restriction, including without limitation
# the rights to use, copy, modify, merge, publish, distribute, sublicense,
# and/or sell copies of the Software, and to permit persons to whom the Software
# is furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included
# in all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
# INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
# IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
# DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
# TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
# OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

"""
labeling.py
"""
import numpy as np
from numpy.random import random_sample

def label(gt_dataset, volume_dim, voxel_dim, labeling_params):
    """
    Labeles the given ground truth dataset according to the config parameters.
    Creates a 3d volume for each fluorophore.

    Args:
        gt_dataset: dict cell_id -> region -> voxels
            the dictionary containing the cell data, splitted by cell_ids and cell regions,
            see load.py for more information.
        volume_dim: (z, x, y) tuplpe
            the dimensions of the ground truth dataset
        voxel_dim: (z, x, y) tuple
            the dimensions of a single voxel in nanometers
        labeling_params: dict
            dictionary containing the labeling parameters for each fluorophore
    Returns:
        labeled_volumes
        labeled_cells

    """
    labeled_volumes = dict()
    labeled_cells = dict()
    for layer in labeling_params:
        fluorophore = labeling_params[layer]['fluorophore']
        del labeling_params[layer]['fluorophore']
        volume, cells  = brainbow(gt_dataset, volume_dim, voxel_dim, **labeling_params[layer])
        if fluorophore in labeled_volumes:
            labeled_volumes[fluorophore] += volume
            labeled_cells[fluorophore] |= cells
        else:
            labeled_volumes[fluorophore] = volume
            labeled_cells[fluorophore] = cells
    return labeled_volumes, labeled_cells


def brainbow(gt_dataset, volume_dim, voxel_dim, region, labeling_density,\
             protein_density, antibody_amp, single_neuron):
    """
    Distributes fluorophores and the corresponding antibodies accross the given cells.
    Follows the brainbow labeling strategy: proteins are distributed
    on the given cell region using a multinomial distribution.
    Then, fluorophore locations are computed by dstirbuting antibodies 
    around the protein locations.

    Args:
    
    Returns:
        labeled_volumes
        labeled_cells

    """
    labeled_cells = {cell_id for cell_id in gt_dataset if random_sample() < labeling_density}
    if single_neuron: 
        labeled_cells = {labeled_cells[0]}
    #Create empty volume
    volume = np.zeros(volume_dim, np.uint32)
    for cell_id in labeled_cells:
        #Get cell data
        reg = list(gt_dataset[cell_id][region])
        if len(reg) == 0: continue
        prob = np.ones(len(reg), np.float64) * 1.0 / len(reg)
        #Compute number of proteins to distribute
        mean_proteins = int(protein_density * len(reg) * np.prod(voxel_dim))
        num_proteins = np.random.poisson(mean_proteins)
        distribution = np.random.multinomial(num_proteins, prob, size=1).squeeze()
        Z, X, Y = np.transpose(reg)
        volume[Z, X, Y] = distribution
        #Add antibodies
        volume = antibodies(volume, voxel_dim, antibody_amp)

    return volume, labeled_cells

def antibodies(volume, voxel_dim, antibody_amp):
    """
    Takes a volume containing binding proteins and distributes antibodies
    accordingly. Adds gaussian noise to mimic the location of random noisy 
    fluorophores in the volume.

    Args:
        volume
        voxel_dim
        antibody_amp
    Returns:
        out
    """
    #Get non zero coordinates
    Z, X, Y = volume.nonzero()
    num_proteins = volume[Z, X, Y]
    #Duplicate voxels
    num_ab = np.round(num_proteins * antibody_amp).astype(np.uint32)
    voxels = np.array((Z, X, Y)).transpose()
    voxels = np.repeat(voxels, num_ab, axis = 0)
    #Add gaussian noise to coordinates
    n = voxels.shape[0]
    ab_std = 0.01 * np.array(voxel_dim)
    noise = np.stack([np.random.normal(0, std, n) for std in ab_std], axis=-1)
    Z, X, Y = np.round(voxels + noise).astype(np.uint32).transpose()
    #Create new volume
    out = np.zeros_like(volume)
    #Ignore out of bound errors
    for i in range(n):
        try:
            out[Z[i], X[i], Y[i]] += 1
        except:
            pass

    return out

