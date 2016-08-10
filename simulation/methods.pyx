# cython: profile=True
import numpy as np
from math import sqrt, log, floor, ceil
import json
import database.methods.exm as exm
import database.methods.datasets as d

cimport numpy as np
cimport database.methods.connectomics as c


#Mishchenko bounds: (542, 813, 93)
#Kasthuri bounds: (10752, 13312, 1849)
#Janelia bounds: (2000, 2000, 6239)

''' Performs expansion of volume in all directions, padding wiht 0s '''
cpdef object perform_expansion(np.ndarray[np.uint32_t, ndim=4] fluo_volume, int expansion_factor):
    cdef np.ndarray[np.uint32_t, ndim=2] number_per_fluor
    (colors, X, Y, Z) = np.nonzero(fluo_volume)
    number_per_fluor = fluo_volume[:, X, Y, Z]
    newX = np.clip(expansion_factor * X + np.random.randint(-expansion_factor, expansion_factor, size = X.size), 0, expansion_factor * fluo_volume.shape[1] - 1)
    newY = np.clip(expansion_factor * Y + np.random.randint(-expansion_factor, expansion_factor, size = X.size), 0, expansion_factor * fluo_volume.shape[2] - 1)
    newZ = np.clip(expansion_factor * Z + np.random.randint(-expansion_factor, expansion_factor, size = X.size), 0, expansion_factor * fluo_volume.shape[3] - 1)
    indices = np.argsort(newZ)
    return (number_per_fluor[:,indices], newX[indices], newY[indices], newZ[indices])

''' Performs expansion of volume in all directions, padding wiht 0s '''
cpdef object perform_expansion_gt(np.ndarray[np.uint32_t, ndim=3] fluo_gt_volume, int expansion_factor):
    (X, Y, Z) = np.nonzero(fluo_gt_volume)
    values = fluo_gt_volume[X, Y, Z]
    newX = expansion_factor * X
    newY = expansion_factor * Y
    newZ = expansion_factor * Z
    indices = np.argsort(newZ)
    
    return (values[indices], newX[indices], newY[indices], newZ[indices])

cpdef object convert_params_to_string(object labeling_unit, object optical_unit):
    labeling_params = labeling_unit.get_param_dict()
    optics_params = optical_unit.get_param_dict()
    parameters = {'labeling_params': labeling_params, 'optics_params': optics_params}
    parameter_string = json.dumps(parameters)
    return parameter_string


''' Runs the simulation given a dataset name, bounds and parameters
cpdef run_sim_with_param_suite(object ds, int gt_dataset_id, object bounds, int z_offset,\
 np.ndarray[np.uint32_t, ndim=1] expansion_factor, object labeling_unit, object optical_unit,\
  object parameters_to_vary, object parameter_ranges,  name = None, comments = ""):

    gt_dataset_name = d.get_name(ds, gt_dataset_id)
    if name is None: name = gt_dataset_name + '_' + str(bounds)
    parameter_string = convert_params_to_string(labeling_unit, optical_unit)
    sim_dataset_id = d.create_simulation_dataset(ds, name, expansion_factor, labeling_unit.get_type(), optical_unit.get_type(), parameter_string, gt_dataset_id, bounds, comments)

    dims = [bounds[0], bounds[1], bounds[2], bounds[3], z_offset - 20, z_offset + 20]
    voxel_dims = np.array(d.get_voxel_dimensions_nm(ds, gt_dataset_id), np.int)
    volume_dim = (bounds[1] - bounds[0], bounds[3] - bounds[2], dims[5] - dims[4])

    optical_unit.calculate_parameters(voxel_dims)
    membranes = c.get_voxels_in_volume_by_cell(ds, gt_dataset_id, dims, include_body = False,  include_membrane = True)
    
    cdef np.ndarray[np.uint8_t, ndim=3] fluorophores, fluorophore_volume
    cdef np.ndarray[np.uint8_t, ndim=2] fluo_slice
    cdef int num_channels
    for param in parameters_to_vary['labeling']:
        param_min, param_max, param_step = parameter_ranges[param]
        while param_min <= param_max:
            setattr(labeling_unit, parameter_name, param_min)
            # LABELING
            # Returns a matrix counting fluorophores in each voxel
            num_channels = len(labeling_unit.fluor_types_used)
            fluorophores = labeling_unit.perform_labeling(membranes, voxel_dims)
            param_min += para_step
            for exp in expansion_factor:
                fluorophore_volume = perform_expansion(fluorophores, exp)
                for param2 in parameters_to_vary['optics']:
                    param_min2, param_max2, param_step2 = parameter_ranges[param2]
                    while param_min2 <= param_max2:
                        setattr(optical_unit, parameter_name, param_min)
                        #EXPANSION
                        #DIVIDE UP TO SAVE RAM
                        for channel in range(num_channels):
                            fluo_slice = optical_unit.resolve_slice(fluorophore_volume, labeling_unit.fluor_types_used, channel)
                            exm.add_slice(ds, sim_dataset_id, z_offset, channel, fluo_slice)
'''

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


            




