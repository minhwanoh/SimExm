'''
labeling.pxd

The BrainbowUnit class
'''

import numpy as np
cimport numpy as np
from database.models.dataset import Dataset,  Fluorset

cdef class BrainbowUnit:
    '''
    BrainbowUnit

    Unit performing the labeling simulation.
    Takes a Dataset object as ground truth and generates volumes of fluorophores in space, with shape (fluorophore, X, Y, Z ).
    '''

    #Attributes
    cdef object gt_dataset
    cdef object fluorset # taken straight from software internal database

    cdef object fluo_volumes  # Used to stack results from multiple calls to label_cells
    cdef object parameters # Used to stack parameters from multiple calls to label_cells
    cdef object labeled_cells  # Used to keep track of labeled cells and region for call to get_ground_truth

    #Methods
    cpdef label_cells(self, region_type = *, fluors = *, protein_density = *, labeling_density = *, antibody_amplification_factor = *, fluor_noise = *, membrane_only = *, single_neuron = *)
        '''
        Labels the cells in the volume withe the given parameter. 
        You can call label_cells multiple times and then use get_labeled_volume to retrived the labeled volume.

        region_type (string) : cell region to label. By default 'full'. To use other regions, make sure they were ingested with the ground truth
        fluors ([string]) : list of fluorphore names to use in the labeling. Each fluorophore is given a seperate channel in the output
        protein_density (float) : density of the protein stain, in number of fluorophore per nm^3. 
        labeling_density (float) : percentage of neurons to label
        antibody_amplification_factor (int) : how much to amplify the protein stain
        fluor_noise (float) : percentage of fluorophore noise to add
        membrane_only (boolean) : if True, labels only the membrane of the indicated region
        single_neuron (boolean) : if True, labels a single random neuron in the dataset
        '''
    
    cpdef object perform_labeling(self, object labeling_dict, object fluors, float protein_density, int antibody_amplification_factor)
        '''
        Performs the labeling, using a dictionary mapping cell_ids to their voxel list. Adds a new fluorophore volume to self.fluo_volumes.
        Helper function for label_cells

        labeling_dict (dict) : dictionnary from cell_id to voxels list (numpy Nx3 array, uint32)
        fluors ([string]) : list of fluorphore names to use in the labeling. Each fluorophore is given a seperate channel in the output
        protein_density (float) : density of the protein stain, in number of fluorophore per nm^3. 
        antibody_amplification_factor (int) : how much to amplify the protein stain
        '''
    
    cpdef np.ndarray[np.uint32_t,ndim=4] add_fluorophore_noise(np.ndarray[np.uint32_t,ndim=4] fluo_volume, int poisson_mean)
        '''Adds random number of fluorophores in the volume. In reconstruciton..'''


    cpdef np.ndarray[np.uint32_t,ndim=2] get_infections(int neuron_number, int num_fluorophores)
        '''
        Returns a 2D list indexed by neurons and giving a boolean array of infections for each neuron 

        neuron_numeber (int) : the number of neurons in the set
        num_fluorophores (int) : the number of different fluorophores that may result from a virus infection
        '''
    cpdef np.ndarray[np.uint32_t, ndim=4] get_labeled_volume()
        ''' Returns the labeled volume as a (fluorophore x X x Y x Z) numpy array, uint32 '''


    cpdef np.ndarray[np.uint32_t, ndim=4] get_ground_truth(membrane_only = *)
         ''' 
         Returns the corresponding ground truth labeling as an (X x Y x Z) numpy array, uint32 

         membrane_only (boolean) : if True, ground truth has only membranes. Tip: use False and use imageJ edge detection if you need the membranes too.
         '''

    cpdef add_param(self, object fluors, float protein_density, float labeling_density,\
      antibody_amplification_factor, float fluor_noise, object region_type, int membrane_only, int single_neuron)
        ''' 
        Adds a new set of parameters to the main parameter dictionary. 
        Each fluorophore is marked with a set of parameters.
         If label_cells is called multiple time with the same fluorophore, the parameters are appended.

        Arguments are the same as for label_cells
        ''' 


    cpdef object get_param_dict()
        ''' Returns a dictionary with the different parameters used throughout the labeling simulation '''


    cpdef object get_type()
        ''' Returns the type of the labeling unit, here Brainbow. '''

    


