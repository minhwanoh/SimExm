'''
labeling.pyx

The BrainbowUnit class
'''

import numpy as np
cimport numpy as np
from math import ceil, floor, log, sqrt, pi, exp
from data.models.dataset import Dataset, Fluorset

cdef class BrainbowUnit:
    '''
    BrainbowUnit

    Unit performing the labeling simulation.
    Takes a Dataset object as ground truth and generates volumes of fluorophores in space, with shape (fluorophore, X, Y, Z ).
    '''

    def __cinit__(self, gt_dataset):
        '''
        Init method, sets attributes

        gt_dataset (Dataset) : Dataset object for the labeling unit to query. See data.models.dataset
        '''
        
        self.gt_dataset = gt_dataset 
        self.fluorset = Fluorset() # taken straight from software internal database

        self.fluo_volumes = [] # Used to stack results from multiple calls to label_cells
        self.parameters = [] # Used to stack parameters from multiple calls to label_cells
        self.labeled_cells = {} # Used to keep track of labeled cells and region for call to get_ground_truth

    cpdef label_cells(self, region_type = 'full', fluors = ['ATTO488'], protein_density = 0.8, labeling_density = 0.2,\
     antibody_amplification_factor = 5, fluor_noise = 0, membrane_only = True, single_neuron = False):
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

        if len(fluors) < 1:
            return
        #Prepare data
        self.num_fluorophores += len(fluors)
        labeling_dict = {}
        cells = self.gt_dataset.get_all_cells()

        if single_neuron: 
            cells = [cells[np.random.randint(0, len(cells))]]
            labeling_density = 1

        for cell in cells:
            if np.random.random_sample() <= self.labeling_density: continue
            if not self.labeled_cells.has_key(cell.get_cell_id()):
                self.labeled_cells[cell.get_cell_id()] = {}
            self.labeled_cells[cell.get_cell_id()][region_type] = 1
            labeling_dict[cell.get_cell_id()] = cell.get_all_regions(region_type, membrane_only = membrane_only)

        self.perform_labeling(labeling_dict, fluors, protein_density, labeling_density, antibody_amplification_factor)

    
    cpdef object perform_labeling(self, object labeling_dict, object fluors, float protein_density, int antibody_amplification_factor):
        '''
        Performs the labeling, using a dictionary mapping cell_ids to their voxel list. Adds a new fluorophore volume to self.fluo_volumes.
        Helper function for label_cells

        labeling_dict (dict) : dictionnary from cell_id to voxels list (numpy Nx3 array, uint32)
        fluors ([string]) : list of fluorphore names to use in the labeling. Each fluorophore is given a seperate channel in the output
        protein_density (float) : density of the protein stain, in number of fluorophore per nm^3. 
        antibody_amplification_factor (int) : how much to amplify the protein stain
        '''

        cdef int x, y, z, n, k, num_neurons, num_proteins, num_fluorophores
        cdef np.ndarray[np.uint32_t, ndim=1] neuron_list, X, Y, Z
        cdef np.ndarray[np.uint32_t, ndim=2] locations, infections
        cdef np.ndarray[np.uint32_t, ndim=4] proteins
        cdef np.ndarray[np.float64_t, ndim= 1] probabilities

        #Get neuron list
        neuron_list = np.array(labeling_dict.keys(), np.uint32)
        num_neurons = neuron_list.size
        num_fluorophores = len(fluors)

        (x, y, z) = self.gt_dataset.get_volume_dim()
        proteins = np.zeros((num_fluorophores, x, y, z), np.uint32)

        #Get infections for each neuron
        infections = self.get_infections(num_neurons, num_fluorophores)
        voxel_dim =  self.gt_dataset.get_voxel_dim()

        for n in range(num_neurons): 
            locations = np.array(labeling_dict[neuron_list[n]], np.uint32)
            (X, Y, Z)= np.transpose(locations)
            probabilities = np.ones(X.size, np.float64) * 1.0/<float>X.size
            for k in range(num_fluorophores):
                if infections[n, k] > 0:
                    num_proteins = np.random.poisson(<int>(self.protein_density * X.size * np.prod(voxel_dim)))
                    protein_distribution = np.random.multinomial(num_proteins, probabilities, size=1).astype(np.uint32)
                    proteins[k, X, Y, Z] =  protein_distribution * self.antibody_amplification_factor

        self.fluo_volumes.append(proteins)

    
    cpdef np.ndarray[np.uint32_t,ndim=4] add_fluorophore_noise(self, np.ndarray[np.uint32_t,ndim=4] fluo_volume, int poisson_mean):
        '''Adds random number of fluorophores in the volume. In reconstruciton..'''

        cdef int x, y, z, channels
        channels, x, y, z = fluo_volume.shape[0], fluo_volume.shape[1], fluo_volume.shape[2], fluo_volume.shape[3]
        cdef np.ndarray[np.uint32_t, ndim=4] poisson_number = np.random.poisson(poisson_mean, size=(channels, x, y, z)).astype(np.uint32)
        cdef np.ndarray[np.uint32_t, ndim=4] poisson_probability = (np.random.rand(channels, x, y, z) < self.fluor_noise).astype(np.uint32)
        return np.add(fluo_volume, np.multiply(poisson_number, poisson_probability) * self.antibody_amplification_factor)

    cpdef np.ndarray[np.uint32_t,ndim=2] get_infections(self, int neuron_number, int num_fluorophores):
        '''
        Returns a 2D list indexed by neurons and giving a boolean array of infections for each neuron 

        neuron_numeber (int) : the number of neurons in the set
        num_fluorophores (int) : the number of different fluorophores that may result from a virus infection
        '''

        cdef int i, j
        cdef np.ndarray[np.uint32_t,ndim=2] infect = np.zeros((neuron_number, num_fluorophores), np.uint32)
        for i in range(neuron_number):
                while np.max(infect[i,:]) == 0:
                    for j in range(num_fluorophores):
                        infect[i,j] = np.random.random_sample() < 0.5
        return infect

    cpdef cdef np.ndarray[np.uint32_t, ndim=4] get_labeled_volume(self):
        ''' Returns the labeled volume as a (fluorophore x X x Y x Z) numpy array, uint32 '''

        return np.concatenate(self.fluo_volumes, 0)


    cpdef cdef np.ndarray[np.uint32_t, ndim=4] get_ground_truth(self, membrane_only = False):
         ''' 
         Returns the corresponding ground truth labeling as an (X x Y x Z) numpy array, uint32 

         membrane_only (boolean) : if True, ground truth has only membranes. Tip: use False and use imageJ edge detection if you need the membranes too.
         '''

        cdef int x, y, z
        (x, y, z) = self.gt_dataset.get_volume_dim()
        cdef np.ndarray[np.uint32_t, ndim=4] ground_truth = np.zeros((x, y, z), np.uint32)
        cells = self.gt_dataset.get_all_cells():

        for cell in cells:
            if self.labeled_cells.has_key(cell_id):
                for region in self.labeled_cells[cell_id].keys():
                    locations = cell.get_all_regions_of_type(region, membrane_only = membrane_only)
                    (X, Y, Z) = np.transpose(locations)
                    ground_truth[X, Y, Z] = cell_id

        return ground_truth

    cpdef add_param(self, int num_fluorophores protein_density, labeling_density, antibody_amplification_factor, fluors, fluor_noise):
        ''' Returns a dicitonary containing the parameters used by the labeling simulation ''' 

        params = {}
        params['protein_density'] = self.protein_density 
        params['labeling_density']  = self.labeling_density 
        params['antibody_amplification_factor'] = self.antibody_amplification_factor 
        params['fluor_types_used'] = self.fluor_types_used
        params['num_fluorophores'] = self.num_fluorophores
        params['fluor_noise'] = self.fluor_noise
        return params

    cpdef get_param_dict(self)

    cpdef get_type(self):
        return "Brainbow"
    


