'''
lib.py

Library for the database module
'''

import numpy as np
from models.dataset import Dataset
from models.cell import Cell, CellRegion
from scipy.ndarray import imread

def load_to_array(image_path, bounds):
	'''
	Loads an image_sequence laocated at image_path into a numpy volume  of shape:
	(bounds[5] - bounds[4], bounds[1] - bounds[0], bound[3] - bounds[2]), uint32
	Images should be named in an ordered sequence and the folder must not contain any other files)

	image_path (string) : the path from which to load the images
	bounds (tuple (x1, x2, y1, y2, z1, z2)) : the bounds to load from the images
	'''

	if image_path[-1] != '/': image_path += '/'
	try:
		os.remove(image_path + '.DS_Store')
	except:
		pass

	im_data = np.zeros((bounds[5] - bounds[4], bounds[1] - bounds[0], bound[3] - bounds[2]), np.uint32)
	images = os.listdir(image_path)
	for i in range(bounds[5] - bounds[4]):
		im_data[i,:,:] = imread(image_path + image, mode = 'L')

	return im_data

	
def  get_composite(cell_dict, dim, max_number_cell = None):
	''' 
	Create composite from list of neurons in a given volume. Can be useful to visualize ground truth data.
	Returns a volume of shape (z, x, y) where x, y, z are taken from the dim argument.

	cell_dict (dict) : dict from cell_id to voxel locations
	dim (tuple) : the (x, y, z) dimension of the volume that the cells reside in
	max_number_cell (int) : if set, limits the number of cell displayed to this number
	'''

    cell_ids = cell_dict.keys()
    (x, y, z) = dim
    volume = np.zeros((z, x, y, 3), np.uint8)
    bound = min(len(cell_ids), max_number_cell) if max_number_cell is not None else len(cell_ids)
    for i in range(bound):
        cell_id = cell_ids[i]
        (Xsm, Ysm, Zsm) = np.split(cell_dict[cell_id], 3, axis=1)
        colorm = [np.random.randint(0, 256) for i in xrange(3)]
        volume[Zsm, Xsm, Ysm, :] = colorm
    return volume



