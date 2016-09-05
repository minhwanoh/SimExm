path_to_SimExm = '/home/jeremy/SimExm'
import sys; sys.path.insert(0, path_to_SimExm)

import os
from src.database.lib import load_to_array, get_composite
from src.database.models.dataset import Dataset
from src.database.models.sim import SimStack, SimParams
from src.simulation.models.expansion import ExpansionUnit
from src.simulation.models.labeling import BrainbowUnit
from src.simulation.models.optics import ConfocalUnit
from PIL import Image
import numpy as np
from src.simulation._psf import gaussian_sigma



def test_data_input(path):

	gt_dataset = Dataset()
	offset = (0, 0, 0)
	bounds_required = [1, 2000, 2000]
	im_data = load_to_array(path, offset, bounds_required)

	gt_dataset.load_from_image_stack((8,8,8), im_data)

	print gt_dataset.get_all_cell_ids()

	cell =  gt_dataset.get_all_cells()['133']

	print cell.get_region_types()
	print cell.get_all_region_ids('full')
	print cell.get_region('full', 1)
	print cell.get_all_regions_as_array('full')

	print gt_dataset.get_voxel_dim()
	print gt_dataset.get_volume_dim()


def view_gt(path, slice_id):

	gt_dataset = Dataset()
	offset = (slice_id, 0, 0)
	bounds_required = [1, 2000, 2000]

	im_data = load_to_array(path, offset, bounds_required)
	gt_dataset.load_from_image_stack((8,8,8), im_data)
	cells = gt_dataset.get_all_cells()

	label_dict = {cell_id : cells[cell_id].get_full_cell(membrane_only=True) for cell_id in cells.keys()}
	im_array = np.squeeze(get_composite(label_dict, (1, 2000, 2000)))
	im = Image.fromarray(im_array, 'RGB')
	im.save('/home/jeremy/test.png')

def get_psf():
    sigma = gaussian_sigma(0.488, 0.488, 1.15, 1.33, 0.1, widefield = False, paraxial=False)
    
    print sigma

#test_data_input("/home/jeremy/connectomics_data/ground_truth_original_format/Janelia/images")
#view_gt("/home/jeremy/connectomics_data/ground_truth_original_format/Janelia/images", 0)

get_psf()