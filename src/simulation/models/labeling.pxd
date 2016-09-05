'''
labeling.pxd

The BrainbowUnit class
'''

import numpy as np
cimport numpy as np

cdef class BrainbowUnit:
    '''
    BrainbowUnit

    Unit performing the labeling simulation.
    Takes a Dataset object as ground truth and generates volumes of fluorophores in space, with shape (fluorophore, X, Y, Z ).
    '''

    #Attributes
    cdef object gt_dataset

    cdef object fluo_volumes  # Used to stack results from multiple calls to label_cells
    cdef object parameters # Used to stack parameters from multiple calls to label_cells
    cdef object labeled_cells  # Used to keep track of labeled cells and region for call to get_ground_truth

    #Methods
    cpdef label_cells(self, object region_type=*, object fluors=*, float  protein_density=*, float labeling_density=*,\
     int antibody_amplification_factor=*, float fluor_noise=*, membrane_only=*, single_neuron=*)
    
    cpdef object perform_labeling(self, object labeling_dict, object fluors, float protein_density, int antibody_amplification_factor)
    
    cpdef np.ndarray[np.uint32_t,ndim=4] add_fluorophore_noise(self, np.ndarray[np.uint32_t,ndim=4] fluo_volume, int poisson_mean)

    cpdef np.ndarray[np.uint32_t,ndim=2] get_infections(self, int neuron_number, int num_fluorophores)

    cpdef np.ndarray[np.uint32_t, ndim=4] get_labeled_volume(self)

    cpdef np.ndarray[np.uint32_t, ndim=4] get_ground_truth(self, membrane_only = *)

    cpdef add_param(self, object fluors, float protein_density, float labeling_density,\
      antibody_amplification_factor, float fluor_noise, object region_type, int membrane_only, int single_neuron)

    cpdef object get_parameters(self)

    cpdef object get_fluors_used(self)

    cpdef object get_type(self)
    


