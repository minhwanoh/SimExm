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

    expansion_factor (int) : the expansion factor to use, currently requires an integer
    '''

    #Attributes
    cdef int expansion_factor

    #Methods
    cpdef object perform_expansion(np.ndarray[np.uint32_t, ndim=4] fluo_volume)
        ''' 
        Performs expansion of volume in all directions, padding wiht 0s 

        fluo_volume (fluorophore x X x Y x Z numpy array, uint32) : the volume to expand
        '''

    cpdef object perform_expansion_gt(np.ndarray[np.uint32_t, ndim=3] fluo_gt_volume)
        ''' 
        Performs expansion of the groudn truth volume in all directions, padding wiht 0s 

        fluo_volume_gt (X x Y x Z numpy array, uint32) : the ground truth volume to expand
        '''

    cpdef object get_param_dict(self)
        ''' Returns dictionary containing the expansion parameters '''