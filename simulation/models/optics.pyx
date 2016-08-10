# cython: profile=False
import numpy as np
import cython
cimport cython
cimport numpy as np
import scipy.ndimage
import scipy.misc
import scipy.signal
from math import ceil, floor, log, sqrt, pi, exp, asin
import database.methods.fluorophore as f
from cython.parallel cimport prange, parallel
from visualization.methods import make_gif


##################### OPTICS #####################

cdef class Confocal:
    def __cinit__(self, ds, num_channels =3,  laser_wavelengths = [488, 565, 660], laser_power  = [50, 50, 50], \
            laser_percentage   = [0.25, 0.25, 0.25], objective_back_aperture  = 1.0, baseline = [30,30,30],\
        std_dev_baseline = [50, 50, 50], filters = [[450, 550], [550, 625], [625, 950]], exposure_time  = 0.1 ,\
        numerical_aperture  = 1.15, objective_efficiency  = 0.8, detector_efficiency  = 0.6,\
        focal_plane_depth  = 500, objective_factor  = 40.0, pixel_size = 6500, precision  = 1.0):
        self.ds = ds
        self.num_channels = num_channels
        self.laser_wavelengths = laser_wavelengths # [nm]
        self.laser_power  = laser_power # [milliWatts]
        self.laser_percentage   = laser_percentage# 1x3 fraction
        self.objective_back_aperture  = objective_back_aperture # [cm]
        self.baseline = baseline# number of photons addew   qqd to everything
        self.std_dev_baseline = baseline # standard deviation of baseline
        self.filters = filters# [nm]
        self.exposure_time  = exposure_time # [seconds]
        self.numerical_aperture  = numerical_aperture
        self.objective_efficiency  = objective_efficiency
        self.detector_efficiency   = detector_efficiency
        self.focal_plane_depth  = focal_plane_depth # [nm] thickness of each focus plane
        self.objective_factor  = objective_factor # magnification factor of objective lens
        self.pixel_size  = pixel_size # [nm] the width of each pixel
        self.precision  = precision # internal parameter used in blurring)
        #assert((num_channels == length(laser_wavelengths)) and (num_channels == length(laser_intensities)) and (num_channels == filters_size(1)))

    ''' Calculates remaining parameters '''
    cpdef calculate_parameters(self, object voxel_dims, int num_slices):
        self.laser_radius = self.objective_back_aperture / self.objective_factor
        self.laser_intensities = [(self.laser_power[i] * self.laser_percentage[i]) * 1.0 / (pi * (self.laser_radius) ** 2) for i in range(len(self.laser_power))]# [Watts / (m ^ 2)]
        self.scale_factor_xy  = voxel_dims[0] * self.objective_factor * self.precision / self.pixel_size
        self.z_offset_step = <int>ceil(<float>self.focal_plane_depth / <float>voxel_dims[2])
        self.std_dev_xy = [(self.scale_factor_xy * self.laser_wavelengths[i]) / (4.44 * self.numerical_aperture * voxel_dims[0]) for i in range(len(self.laser_wavelengths))]
        self.std_dev_z  = [((1.037 * self.laser_wavelengths[i])  * self.precision) / ((self.numerical_aperture ** 2) *  voxel_dims[2])  for i in range(len(self.laser_wavelengths))]
        cdef int slices_required = <int> ceil(8 * max(self.std_dev_z) + <float>(num_slices * self.focal_plane_depth) / <float> voxel_dims[2])
        return slices_required

    ''' Main optics function, resolves fluorophore volume given the microscope parameters ''' 
    cpdef np.ndarray[np.uint8_t, ndim=3] resolve_volume(self, object volume, object volume_dims, object all_fluor_types, int channel):

        cdef int x, y, z, offset, i, start_index, end_index
        cdef np.ndarray[np.uint32_t, ndim = 1] mean_photons_per_fluor, photons
        cdef np.ndarray[np.float64_t, ndim = 3] non_normalized_volume
        cdef np.ndarray[np.uint8_t, ndim=3] normalized_volume
        
        (num_fluors_per_channel, X, Y, Z) = volume

        #Perform magnification
        x = ceil(ceil(volume_dims[0] * self.scale_factor_xy) * (1.0 / self.precision))
        y = ceil(ceil(volume_dims[1] * self.scale_factor_xy) * (1.0 / self.precision))
        X = np.floor(np.floor(self.scale_factor_xy * X) * (1.0 / self.precision)).astype(np.int64)
        Y = np.floor(np.floor(self.scale_factor_xy * Y) * (1.0 / self.precision)).astype(np.int64)

        #Get initial photon counts
        mean_photons_per_fluor = self.get_mean_number_photons_per_fluor(all_fluor_types, channel)
        photons = self.get_photon_count(np.array(num_fluors_per_channel, np.uint32), mean_photons_per_fluor)

        #Resolve slices 1 by 1
        i = 0
        start_index = int(ceil(4 * max(self.std_dev_z)))
        end_index = int(floor(volume_dims[2] - 4 * max(self.std_dev_z)))
        non_normalized_volume = np.zeros(((end_index - start_index) / self.z_offset_step + 1, x, y), np.float64)

        for offset in range(start_index, end_index, self.z_offset_step):
            non_normalized_volume[i,:,:] = self.resolve_slice(X, Y, Z, photons, (x,y), offset, channel)
            i = i + 1

        #Normalize
        normalized_volume = self.normalize(non_normalized_volume)

        return normalized_volume

    ''' Creates the corresponding ground truth volume '''
    cpdef np.ndarray[np.uint32_t, ndim = 3] resolve_ground_truth(self, object volume_gt, object volume_dims):

        cdef int z_low, z_high, x, y, z, i
        (values, X, Y, Z) = volume_gt
        #Perform magnification
        x = ceil(ceil(volume_dims[0] * self.scale_factor_xy) * (1.0 / self.precision))
        y = ceil(ceil(volume_dims[1] * self.scale_factor_xy) * (1.0 / self.precision))
        X = np.floor(np.floor(self.scale_factor_xy * X) * (1.0 / self.precision)).astype(np.int64)
        Y = np.floor(np.floor(self.scale_factor_xy * Y) * (1.0 / self.precision)).astype(np.int64)
        #Create ground truth volume

        i = 0
        start_index = int(ceil(4 * max(self.std_dev_z)))
        end_index = int(floor(volume_dims[2] - 4 * max(self.std_dev_z)))
        cdef np.ndarray[np.uint32_t, ndim = 3] ground_truth = np.zeros(((end_index - start_index) / self.z_offset_step + 1, x, y), np.uint32)
        for offset in range(start_index, end_index, self.z_offset_step):
            z_low = np.searchsorted(Z, offset, side = 'left')
            z_high = np.searchsorted(Z, Z[z_low], side = 'right') - 1
            ground_truth[i, X[z_low: z_high], Y[z_low: z_high]] = values[z_low: z_high]
            i += 1

        return ground_truth


    ''' Resolves the slice optically, given the number and type of fluoropores per voxel'''
    cpdef np.ndarray[np.float64_t, ndim=2] resolve_slice(self, object X, object Y, object Z, object photons, object dims, int offset, int channel):

        cdef int z_low, z_high, x , y
        cdef np.ndarray[np.float64_t, ndim = 2] non_normalized

        (x, y) = (dims[0], dims[1])
        z_low = np.searchsorted(Z, offset - 4 * self.std_dev_z[channel],  side = 'left')
        z_high = np.searchsorted(Z, offset + 4 * self.std_dev_z[channel], side = 'right')
        non_normalized = self.project_photons(X[z_low: z_high], Y[z_low: z_high], Z[z_low: z_high],\
         photons[z_low: z_high], offset, channel, (x, y))

        #Resolve
        non_normalized = np.add(non_normalized, self.get_baseline_image(non_normalized.shape[0], non_normalized.shape[1], channel))

        return non_normalized

    ''' Calculates the number of photons per fluor for a given channel and returns a 1D array of length len(fluors) '''
    cdef np.ndarray[np.uint32_t, ndim=1] get_mean_number_photons_per_fluor(self, object fluors, int channel):
        cdef int i, emitted_photons, detected_photons
        cdef int num_channels = len(fluors)
        cdef np.ndarray[np.uint32_t, ndim=1] mean_detected_photons = np.zeros(num_channels, np.uint32)

        for i in range(num_channels):
            fluor = fluors[i]
            emitted_photons = self.get_emitted_photons(fluor, channel)#Number of emitted photons
            mean_detected_photons[i] = self.get_detected_photons(fluor, emitted_photons, channel)#Number of detected photons

        return mean_detected_photons

    ''' Returns the number of emitted photons given a fluor and the channel '''
    cdef int get_emitted_photons(self, object fluor, int channel):
        cdef float quantum_yield = f.get_quantum_yield(self.ds, fluor)
        cdef float extinction_coefficient = f.get_extinction_coefficient(self.ds, fluor)
        cdef float excitation = f.find_excitation(self.ds, fluor, self.laser_wavelengths[channel])
        cdef float CONSTANT = 0.119626566 # in m^3 * kg * s^{-1} Avogadro's number * Planck's constant * speed of light
        cdef float photons = excitation * quantum_yield * (extinction_coefficient * 100) * self.exposure_time\
         * self.laser_intensities[channel] * (self.laser_wavelengths[channel] * 1e-9) / (1000 * CONSTANT)

        return <int>(ceil(photons))

    ''' Returns the number of detected photons given the fluor and the number of emitted_photons '''
    cdef int get_detected_photons(self, object fluor, int emitted_photons, int channel):
        cdef float wavelength_min = self.filters[channel][0]
        cdef float wavelength_max = self.filters[channel][1]
        cdef float emission = f.find_emission(self.ds, fluor, wavelength_min, wavelength_max)
        cdef float detected_photons = emitted_photons * emission * self.objective_efficiency * self.detector_efficiency
        return <int>(ceil(detected_photons))

    ''' Takes a list of fluorophore locations and the number per type of fluor, and returns the photon count for each of these locations'''
    cdef np.ndarray[np.uint32_t, ndim = 1] get_photon_count(self, np.ndarray[np.uint32_t, ndim=2] num_fluors_per_channel,\
     np.ndarray[np.uint32_t, ndim=1] mean_detected_photons):

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

    ''' Get baseline image for initial photon count '''
    cdef np.ndarray[np.float64_t, ndim=2] get_baseline_image(self, int x, int y, int channel):
        cdef np.ndarray[np.float64_t, ndim=2] psf = self.get_2dgaussian_kernel(x, y, <float> x / 2.0, <float> y / 2.0) 
        return np.multiply(self.std_dev_baseline[channel] * np.random.uniform(0.2, 0.8, size=(x, y)) + self.baseline[channel], psf)

    ''' Creates point spread image with maximum magnitude 1 '''
    cdef np.ndarray[np.float64_t, ndim=2] get_2dgaussian_kernel(self, int x, int y, float sigma_x, float sigma_y):
        cdef float x_0 = <float>x / 2.0
        cdef float y_0 = <float>y / 2.0
        cdef np.ndarray[np.float64_t, ndim=3] indices = np.indices((x, y), dtype=np.float64)
        cdef np.ndarray[np.float64_t, ndim=2] gaussian = np.exp(-((indices[0] - x_0)**2 / (2* sigma_x ** 2) + (indices[1] - y_0)**2 / (2* sigma_y ** 2)))
        return gaussian

    ''' Creates point spread image with maximum magnitude 1'''
    cdef np.ndarray[np.float64_t, ndim=3] get_convolutional_kernel(self, float sigma_x, float sigma_y, float sigma_z, int channel):
        cdef float x_0, y_0, z_0, w_0, z_r, size_z, size_r
        cdef int x, y, z, k
        cdef np.ndarray[np.float64_t, ndim=4] indices
        cdef np.ndarray[np.float64_t, ndim=3] kernel, w_z
        size_r = <float> self.pixel_size / <float>self.objective_factor
        size_z = ceil(<float> self.focal_plane_depth / <float> self.z_offset_step)
        (x, y, z) = (16 * ceil(sigma_x) , 16 * ceil(sigma_y) , 4 * ceil(sigma_z))
        (x_0, y_0, z_0) = (<float>x / 2.0, <float>y / 2.0,  <float>0.0)
        indices = np.indices((z, x, y), dtype=np.float64)
        w_0 =  <float>self.laser_wavelengths[channel] / (<float> pi * self.numerical_aperture)
        z_r = pi * ((w_0)**2) / <float>self.laser_wavelengths[channel]
        w_z = (w_0) * np.sqrt(1 + (size_z * (indices[0] - z_0) / z_r)**2)
        kernel = (w_0) / (w_z) * np.exp(-((size_r * (indices[1] - x_0))**2 + (size_r * (indices[2] - y_0))**2) / (w_z ** 2))
        return kernel

    ''' Projects the given list of photons on the slice '''
    @cython.boundscheck(False)
    cdef np.ndarray[np.float64_t, ndim=2] project_photons(self, long[:] X,\
      long[:] Y, long[:] Z, unsigned int[:] photons, int z_offset, int channel, object dims):
        cdef long i, index, length, x, y, z, index_x, index_y
        cdef np.ndarray[np.float64_t, ndim=3] kernel = self.get_convolutional_kernel(self.std_dev_xy[channel], self.std_dev_xy[channel], self.std_dev_z[channel], channel)
        cdef double[:,:] photon_slice = np.pad(np.zeros(dims, np.float64), (kernel.shape[1] / 2,), mode ='constant', constant_values = 0)  
        length = X.size
        with nogil, parallel(num_threads=1):
            for i in prange(length, schedule='dynamic'):
                z = Z[i] - z_offset
                if z < 0: z = - z
                for x in range(kernel.shape[1]):
                    for y in range(kernel.shape[2]):
                        index_x = x + X[i]
                        index_y = y + Y[i]
                        photon_slice[index_x, index_y] += kernel[z, x, y] * photons[i]
        return np.array(photon_slice[kernel.shape[1] / 2: - kernel.shape[1] / 2, kernel.shape[2] / 2: - kernel.shape[2] / 2], np.float64)

    ''' Normalize array and convert to uint8 '''
    cdef np.ndarray[np.uint8_t, ndim=3] normalize(self, np.ndarray[np.float64_t, ndim=3] non_normalized):
        cdef float the_max = max(<float> np.amax(non_normalized), 1)
        cdef np.ndarray[np.float64_t, ndim=3] normalized = non_normalized * (255.0 / the_max)
        return np.array(np.floor(normalized), np.uint8)

    ''' Returns a dicitonary containing the parameters used by the optical simulation '''
    cpdef object get_param_dict(self):
        params = {}
        params['laser_wavelengths'] = self.laser_wavelengths # [nm]
        params['laser_power'] = self.laser_power # [milliWatts]
        params['laser_percentage']   = self.laser_percentage# 1x3 fraction
        params['objective_back_aperture']  = self.objective_back_aperture # [cm]
        params['baseline'] = self.baseline# number of photons addew  qqd to everything
        params['std_dev_baseline'] = self.std_dev_baseline # standard deviation of baseline
        params['filters'] = self.filters# [nm]
        params['exposure_time'] = self.exposure_time # [seconds]
        params['numerical_aperture']  = self.numerical_aperture
        params['objective_efficiency']  = self.objective_efficiency
        params['detector_efficiency']   = self.detector_efficiency
        params['focal_plane_depth']  = self.focal_plane_depth # [nm] thickness of each focus plane
        params['objective_factor']  = self.objective_factor # magnification factor of objective lens
        params['pixel_size']  = self.pixel_size # [nm] the width of each pixel
        params['precision']  = self.precision # internal parameter used in blurring)
        return params

    cpdef object get_type(self):
        return "Confocal"