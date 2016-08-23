
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
    cdef int focal_plane_depth
    cdef float objective_factor
    cdef int pixel_size


    cpdef calculate_parameters(self, object voxel_dims)
        '''
        Computes additional optics parameters given the size of a voxel, these are used throughout the imaging simulation.
        Gaussian approximation taken from "Gaussian Approximation offluorescence miscorcope point spread function models" by bo Zhang et al.

        voxel_dims ([int] 1 x 3) : size of a voxel in each of its 3 dimensions, in nm
        '''

    cpdef object resolve_volume(self, object volume, object volume_gt, object volume_dims, object fluors)
        '''
        Main optics function, resolves fluorophore volume given the microscope parameters

        volume (fluor x X x Y x Z tuple) : a representation of the fluorophore volume with locations of fluorophores in 3d as a tranposed matrix.
        volume_dims ([int] 1 x 3) : size of the volume in each of its 3 dimensions, in number of voxels
        fluors ([string]) : list of fluorophore corresponding to each channel in the volume. Can be obtained the the labeling parameter dictionary
        '''


    cpdef np.ndarray[np.uint8_t, ndim=3] resolve_channel(self, object volume, object volume_dims, object all_fluor_types, int channel)
        '''
        Resolves the given channel. Itterates over the fluorophores in the volume, computes the emitted photons, and detected photons per pixel.
        Outputs a normalized volume of n slices, where n is the number of slices requested by the user at ground truth loading time.

        volume (tuple(fluors, X, Y, Z)) : tranposed matrix of the fluorophore volume, with expanded coordinates
        volume_dims (int tuple) : size of the volume in voxels per dimension. Since the volume is tranposed above,  helps  figure out how to place the data in the volume
        all_fluor_types : the fluorophores used in labeling (order should correspond to the order of fluors in volume[0])
        channel (int) : the channel to resolve 
        '''
 
    cpdef np.ndarray[np.uint32_t, ndim = 3] resolve_ground_truth(self, object volume_gt, object volume_dims)
        '''
        Generates the ground truth volume corresponding to the unit's parameters.
        The method is called from resolve_volume but may be use on its own by giving it the fluorophore ground truth volume and dimension

        volume (tuple(X, Y, Z)) : tranposed matrix of the ground truth volume, with expanded coordinates
        volume_dims (int tuple) : size of the volume in voxels per dimension. Since the volume is tranposed above, helps figure out how to place the data in the volume
        '''

    cpdef np.ndarray[np.float64_t, ndim=2] resolve_slice(self, object X, object Y, object Z, object photons, object dims, int offset, int channel)
        ''' 
        Resolves a slice optically, given a list of fluorophore locations in the desired bound of 4 x std.

        X ([int]) : list of locations along the x-axis of the photons
        Y ([int]) : list of locations along the y-axis of the photons
        Z ([int]) : list of locations along the z-axis of the photons
        photons ([int]) : number of photon per location
        dims (tuple (x, y)), the size of the slice to create
        offset (int) : the z -offset where to project the fluorohore's emitted photons
        channel (int) : the channel to resolve
        '''

    cdef np.ndarray[np.uint32_t, ndim=1] get_mean_number_photons_per_fluor(self, object fluors, int channel)
        '''
        Calculates the number of photons per fluor for a given channel and returns a 1D array of length len(fluors)

        fluors ([string]) : list of fluorophore names
        channel (int) :  the channel to compute the photon count on
        '''

    cdef int get_emitted_photons(self, object fluor, int channel)
        ''' 
        Returns the number of emitted photons given a fluor and the channel

        fluor (string) : the fluorophore to query 
        channel (int) : calculate the number photons given this channel
        '''

    cdef int get_detected_photons(self, object fluor, int emitted_photons, int channel)
        '''
        Returns the number of detected photons given the fluor and the number of emitted_photons

        fluor (string) : the fluorophore to query 
        emitted_photons (int) : the number of emitted photons
        channel (int) : calculate the number photons given this channel
        '''

    cdef np.ndarray[np.uint32_t, ndim = 1] get_photon_count(self, np.ndarray[np.uint32_t, ndim=2] num_fluors_per_channel,\
     np.ndarray[np.uint32_t, ndim=1] mean_detected_photons)
        '''
        Takes a list of fluorophore locations and the number per type of fluor, and returns the photon count for each of these locations

        num_fluors_per_channel (numpy fluors x voxels) : number of fluors of a certin type per location
        mean_detected photons (numpy uint32 array) : number of photons detected for each fluorophore
        '''

    cdef np.ndarray[np.float64_t, ndim=2] get_baseline_image(self, int x, int y, int channel)
        ''' 
        Get baseline image for initial photon count

        x (int) : width of the baseline image
        y (int) : height of the baseline image
        channel (int) : the channel on which to base to baseline value
        '''

    cdef np.ndarray[np.float64_t, ndim=2] get_2dgaussian_kernel(self, int x, int y, float sigma_x, float sigma_y)
        '''
        Creates point spread image with maximum magnitude 1

        x (int) : width of the kernel 
        y (int) : height of the kernel
        sigma_x (float) : x-std of the kernel
        sigma_y (float) : y-std of the kernel
        '''

    cdef np.ndarray[np.float64_t, ndim=3] get_convolutional_kernel(self, float sigma_x, float sigma_y, float sigma_z, int channel)
        '''
        Creates point spread image with maximum magnitude 1

        Computes the 3d convolution kernel (point spread function of the micriscope). The std_values are used to bound the kernel toa size that 
        allows to reach up to 4 times the std away from the center. Here only half is taken is z for quicker computation.

        sigma_x (float) : x-std of the kernel
        sigma_y (float) : y-std of the kernel
        sigma_z (float) : z-std of the kernel
        channel (int) : the channel to build the kernel based on
        '''

    cdef np.ndarray[np.float64_t, ndim=2] project_photons(self,  long[:] X, long[:] Y, long[:] Z, unsigned int[:] photons, int z_offset, int channel, object dims)
        ''' 
        Projects the given list of photons on the slice. Uses cython memory views, and can be parallelized if the module is compiled with -fopenmp
        Avoid doing a convolution in a sparse volume by doing parallel photon projection
        
        X ([int]) : memory view, list of locations along the x-axis
        Y ([int]) : memory view, list of locations along the y-axis
        Z ([int]) : memory view, list of locations along the z-axis
        photons ([int]) : memory view, list of photon counts
        offset (int) : the z-offset onto which to project the photons (by taking a slice of the psf)
        channel (int) : the channel to work on
        dim (tuple int) : dimension of the image
        '''

    cdef np.ndarray[np.uint8_t, ndim=3] normalize(self, np.ndarray[np.float64_t, ndim=3] non_normalized)
        '''
        Normalize array by dividing by maximum and convert to uint8

        non_normalized (numpy X x Y x Z float) : the non_normalized volume of photons
        '''

    cpdef object get_param_dict(self)
        ''' Returns a dicitonary containing the parameters used by the optical simulation '''


    cpdef object get_type(self)
        ''' Returns the type of the optical unit '''