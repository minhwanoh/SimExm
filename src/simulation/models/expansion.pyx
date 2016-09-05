'''
expansion.pyx

The ExpansionUnit class
'''

import numpy as np
cimport numpy as np


cdef class ExpansionUnit:
    
    '''
    Unit performing the expansion step of the simulation. Works by multiplying the cooridnates by the expansion factor.
    Some variability is attained by distirbuting fluorophores accross the expanded voxels

    '''

    def __cinit__(self, expansion_factor):
        ''' 
        Init method, sets attributes 

        expansion_factor (int) : the expansion factor to use, currently requires an integer
        '''

        self.expansion_factor = expansion_factor 

    cpdef int get_expansion_factor(self):
        ''' Returns the expansion factor of the unit '''

        return self.expansion_factor

    cpdef object expand_volume(self, np.ndarray[np.uint32_t, ndim=4] fluo_volume):
        ''' 
        Performs expansion of volume in all directions, padding wiht 0s 

        fluo_volume (fluorophore x Z x X x Y numpy array, uint32) : the volume to expand
        '''

        cdef np.ndarray[np.uint32_t, ndim=2] number_per_fluor
        (colors, Z, X, Y) = np.nonzero(fluo_volume)
        number_per_fluor = fluo_volume[:, Z, X, Y]
        newX = np.clip(self.expansion_factor * X + np.random.randint(-self.expansion_factor, self.expansion_factor, size = Z.size), 0, self.expansion_factor * fluo_volume.shape[2] - 1)
        newY = np.clip(self.expansion_factor * Y + np.random.randint(-self.expansion_factor, self.expansion_factor, size = Z.size), 0, self.expansion_factor * fluo_volume.shape[3] - 1)
        newZ = np.clip(self.expansion_factor * Z + np.random.randint(-self.expansion_factor, self.expansion_factor, size = Z.size), 0, self.expansion_factor * fluo_volume.shape[1] - 1)
        indices = np.argsort(newZ)

        return (number_per_fluor[:,indices], newZ[indices], newX[indices], newY[indices])

    cpdef object expand_ground_truth(self, np.ndarray[np.uint32_t, ndim=3] fluo_gt_volume):
        ''' 
        Performs expansion of the groudn truth volume in all directions, padding wiht 0s 

        fluo_volume_gt (Z x X x Y numpy array, uint32) : the ground truth volume to expand
        '''

        (Z, X, Y) = np.nonzero(fluo_gt_volume)
        values = fluo_gt_volume[Z, X, Y]
        newX = self.expansion_factor * X
        newY = self.expansion_factor * Y
        newZ = self.expansion_factor * Z
        indices = np.argsort(newZ)
        
        return (values[indices], newZ[indices], newX[indices], newY[indices])

    cpdef object get_parameters(self):
        ''' Returns dictionary containing the expansion parameters '''
        
        return { 'expansion_factor' : self.expansion_factor }