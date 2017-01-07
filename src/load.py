# Provided under BSD license
# Copyright (c) 2017, Jeremy Wohlwend
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without modification,
# are permitted provided that the following conditions are met:
#
#  - Redistributions of source code must retain the above copyright notice,
#    this list of conditions and the following disclaimer.
#  - Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#  - Neither the name of SimExm nor the names of its contributors may be used
#    to endorse or promote products derived from this software without specific
#    prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
# THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
# PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL JEREMY WOHLWEND BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
# EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

"""
load.py

Set of functions to read ground truth data into simulation format. 
The main function is load_gt which loads the ground truth data and returns
the data in simulation format. The other methods are mainly helpers.
"""

#Could consider using a different convolution library
from scipy.signal import fftconvolve
#Not the fastest image loading library but opencv is annoying to install
#Open to suggestions!
from scipy.misc import imread
import numpy as np
import os


def load_gt(image_path, offset, bounds, regions={}):
    """
    Init method. Reads data from image_path, using the given offset and bounds and 
    loads the data into a cell_id->region->voxel dictionary. Computes membranes for the
    main ground truth segmentation and optionally reads out additional region segmentations.

    Args:
        image_path: string
            the path to the directory containing the main ground truth segmentation.
            The images should be in a sequence (not a TIFF stack), 
            and sorted by their file names 
        offset: int (z, x, y) tuple
            the offset from which to load the data
        bounds: int (z, x, y) tuple
            the size of the ground truth data to load. Could be smaller than the actual size.
        regions: dict region (string) -> image_path (string)
            a simple dictionary pointing form addional region annotaiton to load
            using the name of the region as the key and the image path as value

    Returns:
        gt_dataset: dict cell_id (string) -> region (string) -> voxels (list of (z, x, y) tuples)
            the loaded data, in simulation format, with the cell_ids as keys pointing to 
            a sub dict which points from cell regions to lists of voxels in the form 
            of (z, x, y) tuples, where each tuple is a voxel.
    """
    main_data = load_images(image_path, offset, bounds)
    #Compute cytosol and membrane 
    cytosol, membrane = split(main_data)
    loaded_regions = [cytosol, membrane]
    loaded_region_names = ['cytosol', 'membrane']

    #Add addtional regions
    for name in regions:
        path = regions[name]['region_path']
        data = load_images(path, offset, bounds)
        data = data != 0# Only keep binary information for overlap
        loaded_regions.append(data)
        loaded_region_names.append(name)

    #Load to dict
    gt_dataset = load_cells(main_data, loaded_regions, loaded_region_names)
    return gt_dataset

def parse(image_path):
    """
    Parses the given image directory path by
    fixing the path and ignoring hidden files
    Note that the directory should only contain images

    Args:
        image_path: string
            the image directory from which to read
    Returns:
        images: list of strings
            the list of image names in the given directory
    """
    if image_path[-1] != '/': image_path += '/'
    images = sorted(os.listdir(image_path))
    if images[0] == '.directory':
        images = images[1:]
    if images[0] == '.DS_Store':
        images = images[1:]
    return images

def load_images(image_path, offset, bounds):
    """
    Reads the image sequence located at image path into a 3d array

    Args:
        image_path: string
            the path where the image sequence is located
        offset: (z, x y) tuple
            the offset from which to load the data. For instance an offset of (z, x, y)
            starts at the zth image and uses the (x, y) pixel as top left
        bounds: (z, x, y) tuple
            the size of the data to load. Could be smaller than the full size of data,
            should not be larger.
    Returns:
        gt: numpy uint32 3D array
    """
    if image_path[-1] != '/': image_path += '/'
    gt = np.zeros(bounds, np.uint32)
    images = parse(image_path)
    d, w, h = bounds
    z, x, y = offset
    for i in xrange(d):
        im = imread(image_path + images[z + i], mode = 'I')
        gt[i, :, :] = im[x : x + w, y : y + h]
    return gt

def load_cells(main_gt_data, regions, region_names):
    """
    Given the cell segmentation and region annotations,
    loads the ground truth into simulation format.

    Args:
        main_gt_data: numpy 3D uint 32 array
            the main cell segmentation
        regions: list of numpy boolean arrays
            list of volumes indicating a specific cell region.
            By default, contains "cytosol" and "membrane"
        region_names: list of strings
            the corresponding region names, as indexed in regions
    Returns:
        cells: dict cell_id (string) -> region (string) -> voxels (list of (z, x, y) tuples)
            the loaded data, in simulation format, with the cell_ids as keys pointing to 
            a sub dict which points from cell regions to lists of voxels in the form 
            of (z, x, y) tuples, where each tuple is a voxel.
    """
    cells = {}
    cell_ids = np.unique(main_gt_data)[1:]
    for cell_id in cell_ids:
        cells.setdefault(cell_id, {})
        for region, region_name in zip(regions, region_names):
            indices = np.where(np.multiply(main_gt_data, region) == cell_id)
            cells[cell_id][region_name] = np.transpose(indices)
    return cells


def split(gt):
    """
    Takea a ground truth volume, and splits it between membrane and cytosol
    by convolving with an edge kernel

    Args:
        gt: numpy 3d uint32 array
            the 3d volume to convolve with the edge kernel
    Returns:
        cytosol: numpy 3d boolean array
            volume with 1 if the voxel is in the cytosol of some cell, 0 otherwise
        membrane: numpy 3d boolean array
            volume with 1 if the voxel is on the membrane of some cell, 0 otherwise
    """
    edges = get_edges(gt)
    cytosol = edges == 0
    membrane = edges != 0
    return cytosol, membrane

def edge_kernel():
    """
    Returns a 3d kernel for edge detection in an isotropic context.
    A voxel is an edge if any of its neighbours has a different cell_id
    """
    edge_kernel = - 1.0 * np.ones([3, 3, 3], np.float64)
    edge_kernel[1, 1, 1] = 26.0
    return edge_kernel

def get_edges(array):
    """
    Finds the location of edges by convolving with a 3D kernel.
    An edge is any voxel that has a neighbour with a different value.
    
    Args:
        array: numpy 3D array
            the array to convolve
    Returns:
        out: numpy 3D array
            the edge array, with the same shape as the input array
            edges have a non zero value while non edges are set to 0.
    """
    kernel = edge_kernel()
    conv = fftconvolve(array.astype(np.float64), kernel, 'valid')
    conv = np.round(conv).astype(np.uint32)
    conv = np.pad(conv, ((1,1), (1,1), (1,1)), 'constant')
    return conv







