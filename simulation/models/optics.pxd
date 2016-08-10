import numpy as np
cimport numpy as np

cdef class Confocal:
    #Attributes
    cdef public int num_channels
    cdef object ds
    cdef object laser_wavelengths
    cdef public object laser_power
    cdef object laser_percentage
    cdef float objective_back_aperture
    cdef object baseline
    cdef object std_dev_baseline
    cdef object filters
    cdef float exposure_time
    cdef public float numerical_aperture
    cdef public float objective_efficiency
    cdef float detector_efficiency
    cdef float focal_plane_depth
    cdef float objective_factor
    cdef int pixel_size
    cdef float precision
    cdef float laser_radius 
    cdef object laser_intensities 
    cdef float scale_factor_xy 
    cdef int z_offset_step 
    cdef public object std_dev_xy
    cdef public object std_dev_z
    #Methods
    cpdef np.ndarray[np.uint32_t, ndim = 3] resolve_ground_truth(self, object volume_gt, object volume_dims)
    cpdef calculate_parameters(self, object voxel_dims, int num_slices)
    cdef np.ndarray[np.float64_t, ndim=2] get_2dgaussian_kernel(self, int x, int y, float sigma_x, float sigma_y)
    cdef np.ndarray[np.float64_t, ndim=3] get_convolutional_kernel(self, float sigma_x, float sigma_y, float sigma_z, int channel)
    cdef np.ndarray[np.uint32_t, ndim=1] get_mean_number_photons_per_fluor(self, object fluors, int channel)
    cdef int get_emitted_photons(self, object fluor, int channel)
    cdef int get_detected_photons(self, object fluor, int emitted_photons, int channel)
    cdef np.ndarray[np.uint32_t, ndim = 1] get_photon_count(self, np.ndarray[np.uint32_t, ndim=2] num_fluors_per_channel,\
     np.ndarray[np.uint32_t, ndim=1] mean_detected_photons)
    cpdef np.ndarray[np.float64_t, ndim=2] resolve_slice(self, object X, object Y, object Z, object photons, object dims, int offset, int channel)
    cpdef np.ndarray[np.uint8_t, ndim=3] resolve_volume(self, object volume, object volume_dims, object all_fluor_types, int channel)
    cdef np.ndarray[np.float64_t, ndim=2] get_baseline_image(self, int x, int y, int channel)
    cdef np.ndarray[np.float64_t, ndim=2] project_photons(self, long[:] X,\
      long[:] Y, long[:] Z, unsigned int[:] photons, int z_offset, int channel, object dims)
    cdef np.ndarray[np.uint8_t, ndim=3] normalize(self, np.ndarray[np.float64_t, ndim=3] non_normalized)
    cpdef object get_param_dict(self)
    cpdef object get_type(self)