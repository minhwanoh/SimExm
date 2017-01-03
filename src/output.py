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
output.py

Handles storage of simulation outputs.
Simulation stacks can be saved in three possible formats:
    - TIFF stack
    - GIF stack
    - Image sequence (png)
In addition, one can also merge channels into rgb volumes or save each channel
separatly. The same can be done for the ground truth output.
"""

import os
from tifffile import imsave
from optics import scale
import numpy as np

def save_as_tiff(volume, path, name, rgb):
    """
    Saves the given volume at location path/name in TIFF format.
    
    Args:
        volume: 3D array (z, x, y)
            the volume to save
        path: string
            location to use to save the volume
        name: string
            the name of the output volume
        rgb: boolean
            whether to save the volume as an RGB stack or a single channel stack
    """
    dest = path + name + '.tiff'
    if rgb:
        imsave(dest, volume, photometric='rgb')
    else:
        imsave(dest, volume, photometric='minisblack')

def save_as_gif(volume, path, name, rgb):
    """
    Saves the given volume at location path/name in GIF format.
    
    Args:
        volume: 3D array (z, x, y)
            the volume to save
        path: string
            location to use to save the volume
        name: string
            the name of the output volume
        rgb: boolean
            whether to save the volume as an RGB stack or a single channel stack
    """
    dest = path + name + '.gif'
    sequence = [np.squeeze(volume[i]) for i in range(volume.shape[0])]
    writeGif(dest, sequence, duration=0.5)

def save_as_image_sequence(volume, path, name, rgb):
    """
    Saves the given volume at location path/name in a sequence of PNG images,
    with an image for each slice.
    
    Args:
        volume: 3D array (z, x, y)
            the volume to save
        path: string
            location to use to save the volume
        name: string
            the name of the output volume
        rgb: boolean
            whether to save the volume as an RGB stack or a single channel stack
    """
    sequence = [np.squeeze(volume[i]) for i in range(volume.shape[0])]
    digits = np.floor(np.log10(volume.shape[0]))
    for i in range(volume.shape[0]):
        suffix = (digits - np.floor(np.log10(i))) * '0' if i > 0 else digits * '0'
        dest = path + name + '_' + suffix + str(i) + '.png'
        if rgb:
            im = Image.fromarray(np.squeeze(volume[i]), 'RGB')
        else:
            im = Image.fromarray(np.squeeze(volume[i]), 'L')
        im.save(dest)

#The three possible saving methods
SAVE_FUNCTION = {'tiff': save_as_tiff,\
                 'gif': save_as_gif,\
                 'image sequence': save_as_image_sequence}

def save(volumes, path, name, sim_channels, format, **kwargs):
    """
    Saves the simulation stack with the given output parameters.

    Args:
        volumes: list of numpy 3D uint8 arrays
            list of volumes to store, one for each channel
        path: string
            the destination path
        name: string
            name of the simulation experiment
        sim_channels: string ('merged' or 'splitted')
            whether to store each channel in a separate stack or 
            create an RGB stack for every sequence of 3 volumes
        format: string ('tiff', 'gif' or 'image sequence')
            the desired output format
    """
    if path[-1] != "/": path.append("/")

    if not os.path.isdir(path + name):
        os.mkdir(path + name)

    dest = path + name + '/simulation'
    if not os.path.isdir(dest): os.mkdir(dest)
    if sim_channels == 'merged':
        volumes = merge(volumes)
        i = 0
        for vol in volumes:
            sf = SAVE_FUNCTION[format]
            sf(vol, dest + '/', 'channels_{}{}{}'.format(i, i + 1, i + 2), True)
            i += 3
    else:
        i = 0
        for vol in volumes:
            sf = SAVE_FUNCTION[format]
            sf(vol, dest + '/', 'channel_{}'.format(i), False)
            i += 1

def save_gt(cells, labeled_cells, volume_dim, path, name, gt_channels, gt_cells, gt_region, format, optics_params, **kwargs):
    """
    Saves the ground truth stack with the given output parameters
    
    Args:
        labeled_cells
        scaling_factor
        path
        name
        gt_channels
        gt_cells
        gt_region
        format
    """
    if path[-1] != "/": path.append("/")
    if not os.path.isdir(path + name):
        os.mkdir(path + name)

    dest = path + name + '/groundtruth' 
    os.mkdir(dest)
    all_volumes = []
    sf = SAVE_FUNCTION[format]
    for fluorophore in labeled_cells:
        fluorophore_volumes = []
        #Create subdirectory if channels should be splitted
        for cell in labeled_cells[fluorophore]:
            voxels = cells[cell][gt_region]
            volume = np.zeros(volume_dim, np.uint32)
            Z, X, Y = voxels.transpose()
            volume[Z, X, Y] = int(cell)
            fluorophore_volumes.append((volume, cell))
        all_volumes.append(fluorophore_volumes)
        #save ground truth volume
    if gt_channels == 'merged':
        #Merge channels
        volumes = [volume for volume, cell in fluorophore for fluorophore in all_volumes] 
        cells = [cell for volume, cell in fluorophore for fluorophore in all_volumes]
        if gt_cells == 'merged': 
            #Merge cells into a single volume
            volumes = np.sum(volumes, axis=0)
            sf(volumes, dest, 'all_cells', False)
        else:
            #Save each volume separatly
            for volume, cell in zip(volumes, cells):
                sf(volumes, dest, 'cell_' + cell_id, False)
    else:
        for fluorophore, fluorophore_volumes in zip(labeled_cells, all_volumes):
            os.mkdir(dest + '/' + fluorophore)
            volumes = [volume for volume, cell in fluorophore_volumes] 
            cells = [cell for volume, cell in fluorophore_volumes]
            if gt_cells == 'merged':
                #Merge cells
                volume = np.sum(volumes, axis=0)
                sf(volume, dest + '/' + fluorophore, 'all_cells', False)
            else:
                #Save each cell seperatly
                for volume, cell in zip(volumes, cells):
                    sf(cell_volumes, dest + '/' + fluorophore, cell, False)

def merge(volumes):
    """
    Merges the given volumes into multiple 3 channel volumes

    Args:
        volumes: list of numpy 3D arrays
            list of volumes to break up and stack into multiple 3-channels stacks
    Returns:
        out: list of numpy 4D arrays (z, x, y, channel)
            list of rgb volumes
    """
    out = []
    #Add missing channels to round up to % 3
    for i in xrange(3 - len(volumes) % 3):
        empty = np.zeros_like(volumes[0])
        volumes.append(empty)
    #Split every 3 stacks
    for i in xrange(0, len(volumes), 3):
        vol = np.stack(volumes[i:i+3], axis = -1)
        out.append(vol.astype(np.uint8))
    return out