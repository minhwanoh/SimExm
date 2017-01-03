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

Set of methods to handle the labeling of ground truth data.
Implements the brainbow method.
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
        labeled_volumes: dict fluorophore->3d array
            a dict from fluorophore to corresponfing 3d volume
        labeled_cells: dict fluorophore -> list of cell_ids
            list of cells labeled for each fluorophore
    """
    labeled_volumes = dict()
    labeled_cells = dict()
    for layer in labeling_params:
        fluorophore = labeling_params[layer]['fluorophore']
        volume, cells  = brainbow(gt_dataset, volume_dim, voxel_dim, **labeling_params[layer])
        if fluorophore in labeled_volumes:
            labeled_volumes[fluorophore] += volume
            labeled_cells[fluorophore] |= cells
        else:
            labeled_volumes[fluorophore] = volume
            labeled_cells[fluorophore] = cells
    return labeled_volumes, labeled_cells


def brainbow(gt_dataset, volume_dim, voxel_dim, region, labeling_density,\
             protein_density, protein_noise, antibody_amp, single_neuron, **kwargs):
    """
    Distributes fluorophores and the corresponding antibodies accross the given cells.
    Follows the brainbow labeling strategy: proteins are distributed
    on the given cell region using a multinomial distribution.
    Then, fluorophore locations are computed by dstirbuting antibodies 
    around the protein locations.

    Args:  
        gt_dataset: dict cell_id -> region -> voxels
            ground truth dataset in dict format
        volume_dim: (z, x, y) tuple
            dimensions of the ground truth volume
        voxel_dim: (z, x, y) tuple
            dimensions of a ground truth voxel
        region: list of (z, x, y) tuples
            list of voxels
        labeling_density: float64
            the proportion of cells to label
       protein_density: float
            the density of protein labeling, in protein per nm^3
       protein_noise: float
            the amount of protein noise to include.
            Determines what proportion of proteins flies away
            from the labeled region
       antibody_amp: float
            the factor by which to amplify the protein density
       single_neuron: boolean
            if True, only a single cell is labeled
    
    Returns:
        labeled_volumes: numpy uint32 3D array
            the labeled volume
        labeled_cells: set
            the set of cell_ids, indicating which cells were labeled
    """
    labeled_cells = {cell_id for cell_id in gt_dataset if random_sample() < labeling_density}
    if single_neuron: 
        labeled_cells = {labeled_cells[0]}
    #Create empty volume
    volume = np.zeros(volume_dim, np.uint32)
    for cell_id in labeled_cells:
        #Get cell data
        reg = gt_dataset[cell_id][region]
        if len(reg) == 0: continue
        prob = np.ones(len(reg), np.float64) * 1.0 / len(reg)
        #Compute number of proteins to distribute
        mean_proteins = int(protein_density * len(reg) * np.prod(voxel_dim))
        num_proteins = np.random.poisson(mean_proteins)
        distribution = np.random.multinomial(num_proteins, prob, size=1).squeeze()
        voxels = noise(reg, volume_dim, voxel_dim, protein_noise)
        Z, X, Y = voxels.transpose()
        for i in range(len(Z)):
            volume[Z[i], X[i], Y[i]] += np.round(distribution[i] * antibody_amp).astype(np.uint32)

    return volume, labeled_cells

def noise(voxels, volume_dim, voxel_dim, protein_noise):
    """
    Adds gaussian noise to a random subset of the given voxels.
    Clips values outside of the volume dimensions.

    Args:
        voxels: numpy 2D array (n x 3)
            list of (z, x, y) tuple representing a voxel
        volume_dim: (z, x, y) tuple
            dimensions of the ground truth volume
        voxel_dim: (z, x, y) tuple
            dimensions of a ground truth voxel
        protein_noise: float
            the amount of protein noise to include.
            Determines what proportion of proteins flies away
            from the labeled region
    Returns:
        voxels: numpy 2D array (n x 3)
            list of voxels including gaussian noise 
    """
    ab_std = 200.0 / np.array(voxel_dim)#200 nm seems to work well as std
    gaussian = np.stack([np.random.normal(0, std, len(voxels)) for std in ab_std], axis=-1)
    #Get random subset
    indices = np.arange(len(voxels))
    np.random.shuffle(indices)
    indices = indices[:int((1 - protein_noise) * len(voxels))]
    #Set others to 0
    gaussian[indices, :] = [0, 0, 0]
    #Add noise
    voxels = np.round(voxels + gaussian).astype(np.uint32)
    #Make sure we don't get voxels which are out of bounds
    voxels = np.clip(voxels, (0,0,0), np.array(volume_dim) - 1)
    return voxels

