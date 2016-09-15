
'''
optics.pxd

The ConfocalUnit class
'''

import cython
import numpy as np
cimport numpy as np

##################### OPTICS #####################

cdef class ConfocalUnit:
    '''
    The ConfocalUnit class. Implements the imaging section of the simulation by replicating a confocal microscope.
    The unit takes a large set of parameters at initialization. It is necessary to call calculate_parameters to get the number of ground truth slices required.

    The method resolve_volume can be called to image the output of the expansion unit. 
    It creates a SimStack object (i.e data.models.sim) which can then be saved to the desired format.


    num_channels (int) : the number of channels to use
    laser_wavelengths ([int]) : list of wavelenghts for the lasers in nm.
    laser_power ([int]) : list of laser power values for each channel in milliWatts
    laser_percentage [float] : value for each channel required (here 1x3 fraction)
    objective_back_aperture (float) : in cm
    baseline ([float]) : baseline number of photons added to the whole image (i.e as the mean of a poisson dist)
    filters ([[int]] channel x 2) : 2 values per channel in nm
    exposure_time (float) : in seconds
    numerical_aperture (float) : the numerical aperture, needed to compute the convolutional kernel
    objective_efficiency (float) : percentage of efficiency of the objecive, 0.0 to 1.0
    detector_efficiency (float) : percentage of efficiency of the detector, 0.0 to 1.0
    focal_plane_depth (int) : the thickness of a z-slice, in nm
    objective_factor (float) : magnification factor of objective lens (20x, 40x, ...)
    pixel_size (int) : the size of a pixel in nm
    '''

    #Attributes
    cdef int num_channels
    cdef object laser_wavelengths
    cdef object laser_power
    cdef object laser_percentage
    cdef float objective_back_aperture
    cdef object baseline
    cdef object filters 
    cdef float exposure_time
    cdef float numerical_aperture
    cdef float objective_efficiency
    cdef float detector_efficiency 
    cdef object laser_radius
    cdef object laser_intensities
    cdef object w_max
    cdef object z_max
    cdef int z_offset_step
    cdef float scale_factor_xy
    cdef int focal_plane_depth
    cdef float objective_factor
    cdef int pixel_size
    cdef float refractory_index
    cdef float pinhole_radius


    cpdef object compute_parameters(self, object voxel_dims, int expansion_factor, object bounds_wanted)

    cpdef object resolve_volume(self, object volume, object volume_dims, object fluors)

    cpdef np.ndarray[np.uint8_t, ndim=3] resolve_channel(self, object volume, object volume_dims, object all_fluor_types, int channel)
 
    cpdef np.ndarray[np.uint32_t, ndim = 3] resolve_ground_truth(self, object volume_gt, object volume_dims)

    cpdef np.ndarray[np.float64_t, ndim=2] resolve_slice(self, object X, object Y, object Z, object photons, object dims, int offset, int channel)

    cdef np.ndarray[np.uint32_t, ndim=1] get_mean_number_photons_per_fluor(self, object fluors, int channel)

    cdef int get_emitted_photons(self, object fluor, int channel)

    cdef int get_detected_photons(self, object fluor, int emitted_photons, int channel)

    cdef np.ndarray[np.uint32_t, ndim = 1] get_photon_count(self, np.ndarray[np.uint32_t, ndim=2] num_fluors_per_channel,\
     np.ndarray[np.uint32_t, ndim=1] mean_detected_photons)

    cdef np.ndarray[np.uint8_t, ndim=3] get_baseline_volume(self, object volume_dim, int channel)

    cdef np.ndarray[np.float64_t, ndim=2] get_2dgaussian_kernel(self, int x, int y, float sigma_x, float sigma_y)

    cdef np.ndarray[np.float64_t, ndim=3] get_convolutional_kernel(self, int channel)

    cdef np.ndarray[np.float64_t, ndim=2] project_photons(self,  long[:] X, long[:] Y, long[:] Z, unsigned int[:] photons, int z_offset, int channel, object dims)

    cdef np.ndarray[np.uint8_t, ndim=3] normalize(self, np.ndarray[np.float64_t, ndim=3] non_normalized)

    cpdef object get_parameters(self)

    cpdef object get_type(self)