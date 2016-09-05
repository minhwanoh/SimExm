'''
lib.py

Library for the database module
'''
import os
import numpy as np
from scipy.misc import imread, imshow

def load_to_array(image_path, offset, bounds):
    '''
    Loads an image_sequence laocated at image_path into a numpy volume  of shape:
    (bounds[0], bounds[1], bounds[2]), uint32
    Images should be named in an ordered sequence and the folder must not contain any other files)

    image_path (string) : the path from which to load the images
    offset (tuple (z, x, y)) : the top left of the output stack
    bounds (tuple (depth, wdith, height)) : the size to crop from the image stack
    '''

    if image_path[-1] != '/': image_path += '/'
    try:
        os.remove(image_path + '.DS_Store')
    except:
        pass

    im_data = np.zeros((bounds[0], bounds[1], bounds[2]), np.uint32)
    images = sorted(os.listdir(image_path))
    im_width, im_height = imread(image_path + images[0], mode = 'I').shape

    assert (offset[1] <= im_width and offset[1] >= 0), "X offset out of bounds"
    assert (offset[2] <= im_height and offset[2] >= 0), "Y offset out of bounds"
    assert (offset[0] <= len(images) and offset[0] >= 0), "Z offset is out of bounds"
    assert (offset[0] + bounds[0] <= len(images)), "Z offset is too high or you are asking for too many slices"
    assert (offset[1] + bounds[1] <= im_width), "Image width is too small to create a simulation output of the requested width"
    assert (offset[2] + bounds[2] <= im_height), "Image height is too small to create a simulation output of the requested height"
    
    for i in range(bounds[0]):
        im_array = imread(image_path + images[offset[0] + i], mode = 'I')
        im_data[i, :, :] = im_array[offset[1]: offset[1] + bounds[1], offset[2]: offset[2] + bounds[2]]

    return im_data

def get_composite(cell_dict, dim, max_number_cell = None):
    ''' 
    Create composite from list of neurons in a given volume. Can be useful to visualize ground truth data.
    Returns a volume of shape (z, x, y) where z, x, y are taken from the dim argument.

    cell_dict (dict) : dict from cell_id to voxel locations
    dim (tuple) : the (z, x, y) dimension of the volume that the cells reside in
    max_number_cell (int) : if set, limits the number of cell displayed to this number
    '''

    cell_ids = cell_dict.keys()
    (z, x, y) = dim
    volume = np.zeros((z, x, y, 3), np.uint8)
    bound = min(len(cell_ids), max_number_cell) if max_number_cell is not None else len(cell_ids)
    for i in range(bound):
        cell_id = cell_ids[i]
        (Zsm, Xsm, Ysm) = np.split(cell_dict[cell_id], 3, axis=1)
        colorm = [np.random.randint(0, 256) for i in xrange(3)]
        volume[Zsm, Xsm, Ysm, :] = colorm

    return volume



