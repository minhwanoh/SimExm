# cython: profile=True
import numpy as np
from math import sqrt, log, floor, ceil
import json
import database.methods.exm as exm
import database.methods.datasets as d

cimport numpy as np
cimport database.methods.connectomics as c


''' Runs the simulation given a dataset name, bounds and parameters'''
cpdef object run_simulation(object ds, int gt_dataset_id, object bounds, int z_offset,\
 int number_of_slices, int expansion_factor, object labeling_unit, object optical_unit, int gt_full):

    cdef np.ndarray[np.uint8_t, ndim=3] fluo_volume, fluorophore_volume
    cdef np.ndarray[np.uint32_t,ndim=4] fluorophores
    cdef int channel, num_channels, i

    voxel_dims = np.array(d.get_voxel_dimensions_nm(ds, gt_dataset_id), np.int)

    cdef int slices_required = optical_unit.calculate_parameters(voxel_dims, number_of_slices)
    slices_required = int(ceil(<float>slices_required / <float>expansion_factor))

    real_bounds = [bounds[0], bounds[1], bounds[2], bounds[3], z_offset - slices_required / 2, z_offset + slices_required / 2]
    dims = [real_bounds[1] - real_bounds[0], real_bounds[3] - real_bounds[2], real_bounds[5] - real_bounds[4]]

    membranes = c.get_voxels_in_volume_by_cell(ds, gt_dataset_id, real_bounds, include_body = False,  include_membrane = True)
    gt_dict = membranes

    if gt_full != 0:
        body = c.get_voxels_in_volume_by_cell(ds, gt_dataset_id, real_bounds, include_body = True,  include_membrane = False)
        gt_dict = body
    # LABELING
    print "Perform labeling..."
    num_channels =  labeling_unit.get_num_channels()
    (fluorophores, fluorophores_gt) = labeling_unit.perform_labeling(membranes, gt_dict, dims, voxel_dims)
    print "Perform expansion..."
    volume = perform_expansion(fluorophores, expansion_factor)
    volume_gt = perform_expansion_gt(fluorophores_gt, expansion_factor)

    #EXPANSION
    #DIVIDE UP TO SAVE RAM
    print "Perform optics..."
    channels = []
    new_dims = [expansion_factor * i for i in dims]
    for channel in range(num_channels):
        fluo_volume = optical_unit.resolve_volume(volume, new_dims, labeling_unit.fluor_types_used, channel)
        channels.append(fluo_volume)

    fluo_volume_gt = optical_unit.resolve_ground_truth(volume_gt, new_dims)

    #fake channel for rgb mode
    channels.append(np.empty_like(fluo_volume)) 

    #Stack channels
    images = []
    images_gt = []
    for i in range(0, fluo_volume.shape[0]):
        images.append(np.stack([np.squeeze(fluorophore_volume[i,:,:]) for fluorophore_volume in channels], axis=-1))
        images_gt.append(fluo_volume_gt[i, :, :])

    return images, images_gt

def convert_params_to_string(object labeling_unit, object optical_unit):
    labeling_params = labeling_unit.get_param_dict()
    optics_params = optical_unit.get_param_dict()
    parameters = {'labeling_params': labeling_params, 'optics_params': optics_params}
    parameter_string = json.dumps(parameters)
    return parameter_string


            




