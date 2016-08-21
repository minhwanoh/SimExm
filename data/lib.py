'''
lib.py

Library for the database module
'''

import numpy as np
from models.dataset import Dataset
from models.cell import Cell, CellRegion


def load_image_stack(db_path, voxel_dim, cell_stack, region_annotations = [], region_types = [], cell_type_dict = None):
	'''
	Loads a ground truth stack into a Dataset object as defined in ./models/dataset.py and returns it.

	db_path (string) : the path to the new db_file
	voxel_dim (tuple) : dimension of a voxel in nm
	cell_stack (numpy X x Y x Z array, uint32) : stack of ground truth data. Each voxel should have a cell_id as value, and 0 if extra-cellular space.
	region_annotations : list of stacks with same shape as cell_stack, where voxels representing certain cell regions are annotated with region ids (for instance a synapse stack)
	region_types : list of strings representing the types of region annotations given. For each stack in region_annotations, region_types stores the corresponding annotated region.
	cell_type_dict (dict {cell_id(int) : type(string)}) : map from cell_id to type (optional)
	'''
	
	dataset = Dataset(db_path, cell_stack.shape, voxel_dim)

	#Find all cell_ids in the volume
	non_zero = cell_stack.non_zero()
	all_cell_ids = cell_stack[non_zero].unique()

	for cell_id in all_cell_ids:
		#Create Cell object
		cell_type = cell_type_dict[cell_id] if cell_type_dict else None
		cell = Cell(cell_id, cell_type)
		#Load main voxels and compute membrane
		load_cell(cell, cell_stack, non_zero)

		#Load other annotations
		for i range(len(region_annotations)):
			#Load region annotations
			load_region(cell, region_annotations[i] region_types[i])

		dataset.add_cell(cell)

	return dataset


def load_cell(cell, cell_stack, non_zero):
	'''
	Helper function to load_image_stack, loads a cell by computing its voxels in three dimension as well as its edges (i.e membrane)

	cell (Cell) : the cell object to which we will add the full and membrane locations
	cell_stack (numpy X x Y x Z array, uint32) : stack of ground truth data. Each voxel should have a cell_id as value, and 0 if extra-cellular space.
	non_zero (numpy X x Y x Z array, bool) : the non_zero locations, for faster loading 
	'''

	#Find cell voxels
	voxel_list = np.where(cell_stack[non_zero] == cell.get_cell_id())
	voxel_list = non_zero[voxel_list]
	(X, Y, Z) = np.tranpose(voxel_list)
	#Hack to compute membranes
	membrane = np.where(cell_stack[X - 1 : X + 1, Y - 1 : Y + 1, Z].unique().size == 1])
	membrane_list = voxel_list[membrane]
	#Set the region 'full', this method is a shortcut (i.e region_id = 1 and region_type = 'full')
	cell.set_full_cell(voxels = voxel_list, membrane = membrane_list)



def load_region(cell, region_stack, region_type):
	'''
	Helper function to load_image_stack. Loads an annotated region. 
	Such region could be synapses, dendrites, axons, mitochondria.. 
	The idea is to use image stacks which have the region_id as value. By superposing the cell_stack,
	we can easily reconstruct which regions belong to what cell and further divide them by region_id, 
    for instance to label a specific synapse and not others.

	cell (Cell) : the cell object to which we will add the region annotation
	region_stack (numpy X x Y x Z array, uint32) : stack of ground truth data. Each voxel should have a region_id as value, and 0 if extra-cellular space.
	region_type (string) : the type of region this data is annotating, for future querying.
	'''

	voxel_list = cell.get_full_cell()
	membrane_list = cell.get_full_cell(membrane_only = True)

	all_region_ids = region_stack[voxel_list].unique()[1:] #remove 0
	for region_id in all_region_ids:
		indices = np.where(region_stack[voxel_list] == region_id)
		membrane_indices = np.where(region_stack[membrane_list] == region_id)
		#create Region object and add to cell
		region = CellRegion(region_id, voxel_list[indices], membrane_list[membrane_indices], region_type)
		cell.add_region(region)

	
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



