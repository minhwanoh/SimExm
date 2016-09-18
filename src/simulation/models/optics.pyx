
'''
optics.pyx

The ConfocalUnit class
'''

import cython
from cython.parallel import parallel, prange
import numpy as np
cimport numpy as np
from math import pi
from src.database.models.dataset import Fluorset
import src.simulation.psf as psf
from src.simulation._psf import gaussian_sigma


cdef class ConfocalUnit:
    '''
    The ConfocalUnit class. Implements the imaging section of the simulation by replicating a confocal microscope.
    The unit takes a large set of parameters at initialization. It is necessary to call calculate_parameters to get the number of ground truth slices required.

    The method resolve_volume can be called to image the output of the expansion unit. 
    It creates a SimStack object (i.e data.models.sim) which can then be saved to the desired format.
    '''

    def __cinit__(self, num_channels = 3,  laser_wavelengths = [488, 565, 660], laser_power  = [50, 50, 50], \
                    laser_percentage   = [0.25, 0.25, 0.25], objective_back_aperture  = 1.0, baseline = [30,30,30],\
                    filters = [[450, 550], [550, 625], [625, 950]], exposure_time  = 0.1 ,\
                    numerical_aperture  = 1.15, objective_efficiency  = 0.8, detector_efficiency  = 0.6,\
                    focal_plane_depth  = 500, objective_factor  = 40.0, pixel_size = 6500, refractory_index = 1.33, pinhole_radius = 0.55):
        '''
        Sets the initial parameters of the confocal unit.

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
        assert (num_channels > 0), "Resolve at least one channel"

        assert (num_channels <= len(laser_wavelengths)),\
         "The number of channels " + str(num_channels) + " is larger than the number of parameters in laser_wavelengths " +str(len(laser_wavelengths))
        assert (num_channels <= len(laser_power)),\
         "The number of channels " + str(num_channels)  + " is larger than the number of parameters in laser_power " +str(len(laser_power))
        assert (num_channels <= len(filters)),\
             "The number of channels " + str(num_channels) + " is larger than the number of parameters in filters " +str(len(filters))
        assert (num_channels <= len(baseline)),\
             "The number of channels " + str(num_channels) + " is larger than the number of parameters in baseline " +str(len(baseline))
        assert (num_channels <= len(laser_percentage)),\
             "The number of channels " + str(num_channels)  + " is larger than the number of parameters in laser_percentage "  +str(len(laser_percentage))

        self.num_channels = num_channels
        self.laser_wavelengths = laser_wavelengths # [nm]
        self.laser_power  = laser_power # [milliWatts]
        self.laser_percentage   = laser_percentage# 1x3 fraction
        self.objective_back_aperture  = objective_back_aperture # [cm]
        self.baseline = baseline # number of photons addew   qqd to everything
        self.filters = filters # [nm]
        self.exposure_time  = exposure_time # [seconds]
        self.numerical_aperture  = numerical_aperture
        self.objective_efficiency  = objective_efficiency
        self.detector_efficiency   = detector_efficiency
        self.focal_plane_depth  = focal_plane_depth # [nm] thickness of each focus plane
        self.objective_factor  = objective_factor # magnification factor of objective lens
        self.pixel_size  = pixel_size # [nm] the width of each pixel
        self.refractory_index = refractory_index #The refractory index of the specimen
        self.pinhole_radius = pinhole_radius #Pinhole radius in um

        self.laser_radius = <float> self.objective_back_aperture /  <float> self.objective_factor
        self.laser_intensities = [(self.laser_power[i] * self.laser_percentage[i]) * 1.0 / (pi * (self.laser_radius) ** 2) for i in range(len(self.laser_power))]# [Watts / (m ^ 2)]
 

    cpdef object compute_parameters(self, object voxel_dims, int expansion_factor, object bounds_wanted):
        '''
        Computes additional optics parameters given the size of a voxel, these are used throughout the imaging simulation.
        The psf model is based on the gaussian beam model. The kernel shape is computed by truncating the psf at 1/e^2 of it's intensity.
        Returns a  tuple (depth, width, height), giving the bounds required from ground truth to achieve the desired output bounds,
        which is function of the confoncal parameters.


        voxel_dims ([int] 1 x 3) : size of a voxel in each of its 3 dimensions, in nm
        expansion_facotr (int) : the expansion factor used in the simulation
        bounds_wanted (tuple (z, x, y) int) : the desired width, height, and depth for the sim output
        
        '''
        self.scale_factor_xy  = <float> (voxel_dims[1] * self.objective_factor)  / <float> self.pixel_size
        self.z_offset_step = np.ceil(<float>self.focal_plane_depth / <float>voxel_dims[0])
        w_0 = [<float>self.laser_wavelengths[i] / (<float> pi * self.numerical_aperture) for i in range(self.num_channels)]
        z_r = [pi * ((w_0[i])**2) / <float>self.laser_wavelengths[i] for i in range(self.num_channels)]
        #e^2 - 1 = 6.389
        self.z_max = [6.389 * z_r[i] / <float> voxel_dims[0] for i in range(self.num_channels)]
        self.w_max = [(<float>self.laser_wavelengths[i] / (pi * self.numerical_aperture))  * np.sqrt(1 + (self.z_max[i] / z_r[i])**2) / voxel_dims[1] for i in range(self.num_channels)]

        cdef int depth_required = np.floor((2 * max(self.z_max) + <float>(bounds_wanted[0] * self.z_offset_step)) /  <float> expansion_factor)
        #TODO : add kernel edge to width, height required
        cdef int width_required = np.ceil(bounds_wanted[1] * 1.0 / (self.scale_factor_xy * expansion_factor))
        cdef int height_required = np.ceil(bounds_wanted[2] * 1.0 / (self.scale_factor_xy * expansion_factor))
        
        return (depth_required, width_required, height_required)

    cpdef object resolve_volume(self, object volume, object volume_dims, object fluors):
        ''' 
        Main optics function, resolves fluorophore volume given the microscope parameters

        volume (fluor x X x Y x Z tuple) : a representation of the fluorophore volume with locations of fluorophores in 3d as a tranposed matrix.
        volume_dims ([int] 1 x 3) : size of the volume in each of its 3 dimensions, in number of voxels
        fluors ([string]) : list of fluorophore corresponding to each channel in the volume. Can be obtained the the labeling parameter dictionary
        ''' 
        cdef np.ndarray[np.uint8_t, ndim=3] fluo_volume
        cdef np.ndarray[np.uint8_t, ndim=4] images

        channels = []
        for channel in range(self.num_channels):
            fluo_volume = self.resolve_channel(volume, volume_dims, fluors, channel)
            channels.append(fluo_volume)

        #Stack channels
        images = np.stack(channels, axis = -1)

        return images


    cpdef np.ndarray[np.uint8_t, ndim=3] resolve_channel(self, object volume, object volume_dims, object all_fluor_types, int channel):
        '''
        Resolves the given channel. Itterates over the fluorophores in the volume, computes the emitted photons, and detected photons per pixel.
        Outputs a normalized volume of n slices, where n is the number of slices requested by the user at ground truth loading time.

        volume (tuple(fluors, Z, X, Y)) : tranposed matrix of the fluorophore volume, with expanded coordinates
        volume_dims (int tuple) : size of the volume in voxels per dimension. Since the volume is tranposed above,  helps  figure out how to place the data in the volume
        all_fluor_types : the fluorophores used in labeling (order should correspond to the order of fluors in volume[0])
        channel (int) : the channel to resolve 
        '''
        cdef int x, y, z, offset, i, start_index, end_index
        cdef np.ndarray[np.uint32_t, ndim = 1] mean_photons_per_fluor, photons
        cdef np.ndarray[np.float64_t, ndim = 3] non_normalized_volume
        cdef np.ndarray[np.uint8_t, ndim=3] normalized_volume
        
        # Load transposed volume
        (num_fluors_per_channel, Z, X, Y) = volume
        num_fluors_per_channel = np.array(num_fluors_per_channel, np.uint32)

        #Perform magnification
        x = np.floor(volume_dims[1] * self.scale_factor_xy)
        y = np.floor(volume_dims[2] * self.scale_factor_xy) 
        X = np.floor(self.scale_factor_xy * X).astype(np.int64)
        Y = np.floor(self.scale_factor_xy * Y).astype(np.int64)

        #Get initial photon counts
        mean_photons_per_fluor = self.get_mean_number_photons_per_fluor(all_fluor_types, channel)
        photons = self.get_photon_count(num_fluors_per_channel, mean_photons_per_fluor)

        i = 0
        start_index = np.ceil(max(self.z_max)) # We start with an offset of 4 times the maximum z-std
        end_index = np.floor(volume_dims[0] - max(self.z_max)) # End with the same offset
        non_normalized_volume = np.zeros(((end_index - start_index) / self.z_offset_step + 1, x, y), np.float64) #Create empty volume

        #Populate volume
        for offset in range(start_index, end_index, self.z_offset_step):
            non_normalized_volume[i,:,:] = self.resolve_slice(Z, X, Y, photons, (x,y), offset, channel)
            i = i + 1

        #Normalize
        normalized_volume = self.normalize(non_normalized_volume)
        #Compiler getting angry if i don't don't unpack manually...
        (vz, vx, vy) = normalized_volume.shape[0], normalized_volume.shape[1], normalized_volume.shape[2]
        normalized_volume = np.add(normalized_volume, self.get_baseline_volume((vz, vx, vy), channel))

        return normalized_volume
 
    cpdef np.ndarray[np.uint32_t, ndim = 3] resolve_ground_truth(self, object volume_gt, object volume_dims):
        '''
        Generates the ground truth volume corresponding to the unit's parameters.
        The method is called from resolve_volume but may be use on its own by giving it the fluorophore ground truth volume and dimension

        volume (tuple(values, Z, X, Y)) : tranposed matrix of the ground truth volume, with expanded coordinates
        volume_dims (int tuple) : size of the volume in voxels per dimension. Since the volume is tranposed above, helps figure out how to place the data in the volume
        '''
        cdef int z_low, z_high, x, y, z, i, start_index, end_index
        (values, Z, X, Y) = volume_gt

        #Perform magnification
        x = np.floor(volume_dims[1] * self.scale_factor_xy)
        y = np.floor(volume_dims[2] * self.scale_factor_xy)
        X = np.floor(self.scale_factor_xy * X).astype(np.int64)
        Y = np.floor(self.scale_factor_xy * Y).astype(np.int64)

        #Create ground truth volume
        i = 0
        start_index = np.ceil(max(self.z_max))
        end_index = np.floor(volume_dims[0] - max(self.z_max))
        cdef np.ndarray[np.uint32_t, ndim = 3] ground_truth = np.zeros(((end_index - start_index) / self.z_offset_step + 1, x, y), np.uint32)

        #Populate volume
        for offset in range(start_index, end_index, self.z_offset_step):
            z_low = np.searchsorted(Z, offset, side = 'left')
            z_high = np.searchsorted(Z, Z[z_low], side = 'right') - 1 if z_low < Z.size else z_low#Corner case
            ground_truth[i, X[z_low: z_high], Y[z_low: z_high]] = values[z_low: z_high]
            i += 1

        return ground_truth

    cpdef np.ndarray[np.float64_t, ndim=2] resolve_slice(self, object Z, object X, object Y, object photons, object dims, int offset, int channel):
        ''' 
        Resolves a slice optically, given a list of fluorophore locations in the desired bound of 4 x std.

        Z ([int]) : list of locations along the z-axis of the photons
        Y ([int]) : list of locations along the y-axis of the photons
        X ([int]) : list of locations along the x-axis of the photons
        photons ([int]) : number of photon per location
        dims (tuple (x, y)), the size of the slice to create
        offset (int) : the z -offset where to project the fluorohore's emitted photons
        channel (int) : the channel to resolve
        '''
        cdef int z_low, z_high, x , y
        cdef np.ndarray[np.float64_t, ndim = 2] non_normalized

        (x, y) = (dims[0], dims[1])
        z_low = np.searchsorted(Z, offset - self.z_max[channel],  side = 'left')
        z_high = np.searchsorted(Z, offset + self.z_max[channel], side = 'right')
        non_normalized = self.project_photons(Z[z_low: z_high], X[z_low: z_high], Y[z_low: z_high], photons[z_low: z_high], offset, channel, (x, y))

        return non_normalized

    cdef np.ndarray[np.uint32_t, ndim=1] get_mean_number_photons_per_fluor(self, object fluors, int channel):
        '''
        Calculates the number of photons per fluor for a given channel and returns a 1D array of length len(fluors)

        fluors ([string]) : list of fluorophore names
        channel (int) :  the channel to compute the photon count on
        '''
        cdef int i, emitted_photons, detected_photons
        cdef int num_channels = len(fluors)
        cdef np.ndarray[np.uint32_t, ndim=1] mean_detected_photons = np.zeros(num_channels, np.uint32)

        for i in range(num_channels):
            fluor = fluors[i]
            emitted_photons = self.get_emitted_photons(fluor, channel)#Number of emitted photons
            mean_detected_photons[i] = self.get_detected_photons(fluor, emitted_photons, channel)#Number of detected photons

        return mean_detected_photons

    cdef int get_emitted_photons(self, object fluor, int channel):
        ''' 
        Returns the number of emitted photons given a fluor and the channel

        fluor (string) : the fluorophore to query 
        channel (int) : calculate the number photons given this channel
        '''
        fluorset = Fluorset()
        f = fluorset.get_fluor(fluor)
        cdef float quantum_yield = f.get_quantum_yield()
        cdef float extinction_coefficient = f.get_extinction_coefficient()
        cdef float excitation = f.find_excitation(self.laser_wavelengths[channel])
        cdef float CONSTANT = 0.119626566 # in m^3 * kg * s^{-1} Avogadro's number * Planck's constant * speed of light
        cdef float photons = excitation * quantum_yield * (extinction_coefficient * 100) * self.exposure_time\
         * self.laser_intensities[channel] * (self.laser_wavelengths[channel] * 1e-9) / (1000 * CONSTANT)

        return np.ceil(photons)

    cdef int get_detected_photons(self, object fluor, int emitted_photons, int channel):
        '''
        Returns the number of detected photons given the fluor and the number of emitted_photons

        fluor (string) : the fluorophore to query 
        emitted_photons (int) : the number of emitted photons
        channel (int) : calculate the number photons given this channel
        '''
        fluorset = Fluorset()
        f = fluorset.get_fluor(fluor)

        cdef float wavelength_min = self.filters[channel][0]
        cdef float wavelength_max = self.filters[channel][1]
        cdef float emission = f.find_emission(wavelength_min, wavelength_max)
        cdef float detected_photons = emitted_photons * emission * self.objective_efficiency * self.detector_efficiency
        return np.ceil(detected_photons)

    cdef np.ndarray[np.uint32_t, ndim = 1] get_photon_count(self, np.ndarray[np.uint32_t, ndim=2] num_fluors_per_channel,\
     np.ndarray[np.uint32_t, ndim=1] mean_detected_photons):
        '''
        Takes a list of fluorophore locations and the number per type of fluor, and returns the photon count for each of these locations

        num_fluors_per_channel (numpy fluors x voxels) : number of fluors of a certin type per location
        mean_detected photons (numpy uint32 array) : number of photons detected for each fluorophore
        '''
        #Initialize with baseline
        cdef int num_voxels = num_fluors_per_channel.shape[1]
        cdef int num_fluors = num_fluors_per_channel.shape[0]

        cdef np.ndarray[np.uint32_t, ndim = 1] mean_photons
        cdef np.ndarray[np.uint32_t, ndim = 1] photons = np.zeros(num_voxels, np.uint32)
        #Add photon count produced by each fluorophore for the given laser
        for fluor in range(num_fluors):
            mean_photons = np.random.poisson(mean_detected_photons[fluor], size=num_voxels).astype(np.uint32) * num_fluors_per_channel[fluor,:]
            photons = np.add(photons, mean_photons)

        return photons

    cdef np.ndarray[np.uint8_t, ndim=3] get_baseline_volume(self, object volume_dim, int channel):
        ''' 
        Creates a baseline volume for initial photon count. 
        Multiplies a random poisson distirbution with mean self.baseline[channel]
        and a gaussian distribution centered at the center of the 2d plane.

        volume_dim (tuple) : dimensions of the volume in pixels
        channel (int) : the channel on which to base to baseline value
        '''
        cdef int x, y, z
        (z, x, y) = volume_dim

        cdef np.ndarray[np.float64_t, ndim=2] psf = self.get_2dgaussian_kernel(x, y, <float> x / 2.0, <float> y / 2.0) 
        cdef np.ndarray[np.uint8_t, ndim=3] baseline_volume = np.zeros(volume_dim, np.uint8)

        for i in range(z):
            baseline_volume[i,:,:] = np.floor(np.multiply(np.random.poisson(self.baseline[channel], size=(x, y)), psf)).astype(np.uint8)
        
        return baseline_volume

    cdef np.ndarray[np.float64_t, ndim=2] get_2dgaussian_kernel(self, int x, int y, float sigma_x, float sigma_y):
        '''
        Creates point spread image with maximum magnitude 1

        x (int) : width of the kernel 
        y (int) : height of the kernel
        sigma_x (float) : x-std of the kernel
        sigma_y (float) : y-std of the kernel
        '''
        cdef float x_0 = <float>x / 2.0
        cdef float y_0 = <float>y / 2.0
        cdef np.ndarray[np.float64_t, ndim=3] indices = np.indices((x, y), dtype=np.float64)
        cdef np.ndarray[np.float64_t, ndim=2] gaussian = np.exp(-((indices[0] - x_0)**2 / (2* sigma_x ** 2) + (indices[1] - y_0)**2 / (2* sigma_y ** 2)))
        return gaussian

    cdef np.ndarray[np.float64_t, ndim=3] get_convolutional_kernel(self,  int channel):
        '''
        Creates point spread image with maximum magnitude 1

        Computes the 3d convolution kernel (point spread function of the micriscope). The std_values are used to bound the kernel toa size that 
        allows to reach up to 4 times the std away from the center. Here only half is taken is z for quicker computation.
        Uses the gaussian beam formula found here: https://en.wikipedia.org/wiki/Gaussian_Beam
        The gaussiam parameters are computed from the external psf library, courtesy of Christophe Golke.

        channel (int) : the channel to build the kernel based on
        '''        
        cdef float x_0, y_0, z_0, w_0, z_r, size_z, size_r
        cdef int x, y, z, k
        cdef np.ndarray[np.float64_t, ndim=4] indices
        cdef np.ndarray[np.float64_t, ndim=3] kernel, w_z

        #Compute size coefficient
        size_r = <float> self.pixel_size / <float>self.objective_factor
        size_z = np.ceil(<float> self.focal_plane_depth / <float> self.z_offset_step)
        
        #Set bounds
        (x, y, z) = (np.ceil(2 * self.w_max[channel]), np.ceil(2 * self.w_max[channel]), np.ceil(self.z_max[channel]))
        (x_0, y_0, z_0) = (<float>x / 2.0, <float>y / 2.0,  <float>0.0)
        indices = np.indices((z, x, y), dtype=np.float64)
        
        #Compute using gaussian beam formula
        w_0 =  <float>self.laser_wavelengths[channel] / (<float> pi * self.numerical_aperture)
        z_r = pi * ((w_0)**2) / <float>self.laser_wavelengths[channel]
        w_z = (w_0) * np.sqrt(1 + (size_z * (indices[0] - z_0) / z_r)**2)
        kernel = (((w_0) / (w_z)) ** 2) * np.exp(- 2 * ((size_r * (indices[1] - x_0))**2 + (size_r * (indices[2] - y_0))**2) / (w_z ** 2))

        #Normalize
        kernel = kernel / np.sum(kernel) 
        
        return kernel

    @cython.boundscheck(False)
    cdef np.ndarray[np.float64_t, ndim=2] project_photons(self,  long[:] Z, long[:] X, long[:] Y, unsigned int[:] photons, int z_offset, int channel, object dims):
        ''' 
        Projects the given list of photons on the slice. Uses cython memory views, and can be parallelized if the module is compiled with -fopenmp
        Avoid doing a convolution in a sparse volume by doing parallel photon projection
        
        Z ([int]) : memory view, list of locations along the z-axis
        X ([int]) : memory view, list of locations along the x-axis
        Y ([int]) : memory view, list of locations along the y-axis
        photons ([int]) : memory view, list of photon counts
        offset (int) : the z-offset onto which to project the photons (by taking a slice of the psf)
        channel (int) : the channel to work on
        dim (tuple int) : dimension of the image
        '''
        cdef long i, index, length, x, y, z, index_x, index_y
        cdef np.ndarray[np.float64_t, ndim=3] kernel = self.get_convolutional_kernel(channel)
        cdef double[:,:] photon_slice = np.zeros(dims, np.float64)
        cdef double[:,:] padded_slice = np.pad(photon_slice, (kernel.shape[1] / 2,), mode ='constant', constant_values = 0)  
        length = Z.size
        with nogil, parallel(num_threads=1):
            for i in prange(length, schedule='dynamic'):
                z = Z[i] - z_offset
                if z < 0: z = - z
                for x in range(kernel.shape[1]):
                    for y in range(kernel.shape[2]):
                        index_x = x + X[i]
                        index_y = y + Y[i]
                        padded_slice[index_x, index_y] += kernel[z, x, y] * photons[i]

        photon_slice[:,:] = padded_slice[kernel.shape[1] / 2 : kernel.shape[1] / 2 + photon_slice.shape[0],\
                            kernel.shape[2] / 2 : kernel.shape[2] / 2 + photon_slice.shape[1]]

        return np.array(photon_slice, np.float64)

    
    cdef np.ndarray[np.uint8_t, ndim=3] normalize(self, np.ndarray[np.float64_t, ndim=3] non_normalized):
        '''
        Normalize array by dividing by maximum and convert to uint8

        non_normalized (numpy Z x X x Y float) : the non_normalized volume of photons
        '''
        cdef float the_max = <float> np.amax(non_normalized)
        cdef np.ndarray[np.float64_t, ndim=3] normalized = non_normalized * (255.0 / the_max)
        return np.array(np.floor(normalized), np.uint8)

    cpdef object get_parameters(self):
        '''Returns a dicitonary containing the parameters used by the optical simulation'''
        params = {}
        params['unit_type'] = self.get_type()
        params['num_channels'] = self.num_channels # number of channles in the optical unit
        params['laser_wavelengths'] = list(self.laser_wavelengths) # [nm]
        params['laser_power'] = list(self.laser_power) # [milliWatts]
        params['laser_percentage']   = list(self.laser_percentage)# 1x3 fraction
        params['objective_back_aperture']  = self.objective_back_aperture # [cm]
        params['baseline'] = list(self.baseline)# number of photons addew  qqd to everything
        params['filters'] = list(self.filters)# [nm]
        params['exposure_time'] = self.exposure_time # [seconds]
        params['numerical_aperture']  = self.numerical_aperture
        params['objective_efficiency']  = self.objective_efficiency
        params['detector_efficiency']   = self.detector_efficiency
        params['focal_plane_depth']  = self.focal_plane_depth # [nm] thickness of each focus plane
        params['objective_factor']  = self.objective_factor # magnification factor of objective lens
        params['pixel_size']  = self.pixel_size # [nm] the width of each pixel

        return params

    cpdef object get_type(self):
        '''Returns the type of the optical unit'''
        return "Confocal"