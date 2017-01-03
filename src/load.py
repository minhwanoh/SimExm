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
load.py

Set of functions to read ground truth data into simulation format. 
The main function is load_gt which loads the groudn truth data and returns
the data in simulation format. The other methods are mainly helpers.
"""

from itertools import product
#Could consider using a different convolution library
from scipy.signal import fftconvolve
#Not the fastest image loading library but opencv is annoying to install
#Open to suggestions!
from scipy.misc import imread
import numpy as np
import os
from matplotlib.pyplot import imshow, show


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
        data = load_images(image_path, offset, bounds)
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
    d, w, h = main_gt_data.shape
    for k, i, j in product(xrange(d), xrange(w), xrange(h)):
        cell_id = main_gt_data[k, i, j]
        if cell_id != 0:
            cells.setdefault(cell_id, {name: [] for name in region_names})
            for region, region_type in zip(regions, region_names): 
                if region[k, i, j] != 0:
                    cells[cell_id][region_type].append((k, i, j))
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
    kernel = edge_kernel()
    edges = np.round(fftconvolve(gt.astype(np.float64), kernel, 'valid'))
    edges = np.pad(edges, ((1,1), (1,1), (1,1)), 'constant')
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





