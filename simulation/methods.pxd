import numpy as np
cimport numpy as np

''' Runs a preview simulation given a dataset name, bounds and parameters'''
cpdef object run_simulation(object ds, int gt_dataset_id, object bounds, int z_offset,\
 int number_of_slices, int expansion_factor, object labeling_unit, object optical_unit, int gt_full)

''' Performs expansion of volume in all directions, padding wiht 0s '''
cpdef object perform_expansion(np.ndarray[np.uint32_t, ndim=4] fluo_volume, int expansion_factor)

cpdef object perform_expansion_gt(np.ndarray[np.uint32_t, ndim=3] fluo_gt_volume, int expansion_factor)

''' Puts the labeling parameters and the otpical parameters in a dictionary for output '''
cpdef object convert_params_to_string(object labeling_unit,object optical_unit)