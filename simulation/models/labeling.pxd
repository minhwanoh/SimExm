import numpy as np
cimport numpy as np

cdef class BrainBowUnit:
    #Attributes
    cdef public float protein_density
    cdef public float labeling_density
    cdef public int antibody_amplification_factor
    cdef public object fluor_types_used
    cdef int num_channels
    cdef public float fluor_noise
    #Methods
    cpdef object perform_labeling(self, object dict_for_labeling, object dict_for_gt, object volume_dims, object voxel_dims)
    cpdef int get_num_channels(self)
    cpdef np.ndarray[np.uint32_t,ndim=4] add_fluorophore_noise(self, np.ndarray[np.uint32_t,ndim=4] fluo_volume, int poisson_mean)
    cpdef np.ndarray[np.uint8_t ,ndim=2] get_infections(self, int neuron_number, int num_channels)
    cpdef get_param_dict(self)
    cpdef get_type(self)