# cython: profile=False
import numpy as np
import cython
cimport numpy as np
from math import ceil, floor, log, sqrt, pi, exp
import database.methods.fluorophore as f

##################### LABELING #####################

cdef class BrainBowUnit:
    def __cinit__(self, protein_density = 0.8, labeling_density = 0.25,\
      fluors = ['ATTO488', 'ATTO550', 'ATTO647N'], antibody_amplification_factor = 5, fluor_noise = 0.00001):
        self.protein_density = protein_density
        self.labeling_density = labeling_density # fraction of cells infected
        self.antibody_amplification_factor = antibody_amplification_factor # to simulate amplification via antibodies
        self.fluor_types_used = fluors # see get_fluor_types.m for a listing of the codes for different types
        self.num_fluorophores = len(fluors)
        self.fluor_noise = fluor_noise

    cpdef int get_num_channels(self):
        return self.num_channels

    ''' Performs the labeling, given a volume containing membrane pixels and using class attributes as parameters. Returns a volume of fluorohore'''
    cpdef object perform_labeling(self, object dict_for_labeling, object dict_for_gt, object volume_dims, object voxel_dims):
        cdef int x, y, z, n, k, num_neurons, num_proteins
        cdef float mean_poisson
        cdef np.ndarray[np.uint32_t, ndim=1] neuron_list, X, Y, Z, gt_X, gt_Y, gt_Z
        cdef np.ndarray[np.uint32_t, ndim=2] locations, gt_locations
        cdef np.ndarray[np.float64_t, ndim= 1] probabilities

        #Prepare volume
        (x, y, z) = volume_dims
        cdef np.ndarray[np.uint32_t, ndim=4] proteins = np.zeros((self.num_channels, x, y, z), np.uint32) # permuting the dimensionality here
        cdef np.ndarray[np.uint32_t, ndim=3] ground_truth = np.zeros((x, y, z), np.uint32)

        if self.num_channels == 0:
            return (proteins, ground_truth)
        #Get neuron list
        neuron_list = np.array(dict_for_labeling.keys(), np.uint32)
        num_neurons = neuron_list.size

        #Get infections for each neuron
        cdef np.ndarray[np.uint8_t, ndim=2] infections = self.get_infections(num_neurons, self.num_channels)

        for n in range(num_neurons): 
            locations = np.array(dict_for_labeling[neuron_list[n]], np.uint32)
            gt_locations = np.array(dict_for_gt[neuron_list[n]], np.uint32) if dict_for_gt.has_key(neuron_list[n]) else np.empty((0, 3), np.uint32)
            (X, Y, Z)= np.transpose(locations)
            (gt_X, gt_Y, gt_Z)= np.transpose(gt_locations)
            probabilities = np.ones(X.size, np.float64) * 1.0/<float>X.size
            for k in range(self.num_channels):
                if infections[n, k] > 0:
                    num_proteins = np.random.poisson(<int>(self.protein_density * X.size))
                    protein_distribution = np.random.multinomial(num_proteins, probabilities, size=1).astype(np.uint32)
                    proteins[k, X, Y, Z] =  protein_distribution * self.antibody_amplification_factor
                    ground_truth[gt_X, gt_Y, gt_Z] = neuron_list[n]

        proteins = self.add_fluorophore_noise(proteins, <int>np.mean(proteins))

        return (proteins, ground_truth)

    '''Adds random number of fluorophores in the volume'''
    cpdef np.ndarray[np.uint32_t,ndim=4] add_fluorophore_noise(self, np.ndarray[np.uint32_t,ndim=4] fluo_volume, int poisson_mean):
        cdef int x, y, z, channels
        channels, x, y, z = fluo_volume.shape[0], fluo_volume.shape[1], fluo_volume.shape[2], fluo_volume.shape[3]
        cdef np.ndarray[np.uint32_t, ndim=4] poisson_number = np.random.poisson(poisson_mean, size=(channels, x, y, z)).astype(np.uint32)
        cdef np.ndarray[np.uint32_t, ndim=4] poisson_probability = (np.random.rand(channels, x, y, z) < self.fluor_noise).astype(np.uint32)
        return np.add(fluo_volume, np.multiply(poisson_number, poisson_probability) * self.antibody_amplification_factor)

    '''Returns a 2D list indexed by neurons and giving a boolean array of infections for each neuron '''
    cpdef np.ndarray[np.uint8_t,ndim=2] get_infections(self, int neuron_number, int num_channels):
        cdef int i, j
        cdef np.ndarray[np.uint8_t,ndim=2] infect = np.zeros((neuron_number, num_channels), np.uint8)
        for i in range(neuron_number):
            if np.random.random_sample() < self.labeling_density:
                while np.max(infect[i,:]) == 0:
                    for j in range(num_channels):
                        infect[i,j] = np.random.random_sample() < 0.5
        return infect

    ''' Returns a dicitonary containing the parameters used by the labeling simulation ''' 
    cpdef get_param_dict(self):
        params = {}
        params['protein_density'] = self.protein_density 
        params['labeling_density']  = self.labeling_density 
        params['antibody_amplification_factor'] = self.antibody_amplification_factor 
        params['fluor_types_used'] = self.fluor_types_used
        params['num_fluorophores'] = self.num_fluorophores
        params['fluor_noise'] = self.fluor_noise
        return params

    cpdef get_type(self):
        return "Brainbow"
    


