'''
dataset.py

The Dataset interface for ground truth data and,
the Fluorset interface to access the Fluorophore objects.

'''
import numpy as np
import h5py
from cell import Cell, CellRegion
import fluors as f


class Dataset:

    '''
    Interface to ingest, manipulate, and store ground truth datasets.
    Labeling units take Datasets as input.

    Uses hdf5 as storage format.
    Use load_from_db and save_to_db as ways to read from or write to a .h5 file

    '''
   
    def __init__(self):
        ''' Init method '''

        self.cells = {}
        self.volume_dim = None
        self.voxel_dim = None

        self.db_path = None
            
    def get_voxel_dim(self):
        ''' Returns the dimensions of a voxel in a 1x3 tuple. Values are in nanometers '''

        return self.voxel_dim

    def get_volume_dim(self):
        ''' Returns the dimensions of the dataset in a 1x3 tuple. Values are in voxel count '''

        return self.volume_dim

    def add_cell(self, cell):
        '''
        Add new Cell object to the dataset

        cell (Cell) : the cell to add. See cell.py
        '''

        #Create cell strucute in the hdf5 format
        self.cells[str(cell.get_cell_id())] = cell

    def get_cell(self, cell_id):
        '''
        Returns the desired Cell object. 

        cell_id (int or string) : the cell_id of the Cell to query.
        '''
        return self.cells[str(cell_id)]

    
    def get_all_cells(self):
        ''' Returns all the cells in the dataset as dictionary from cell_ids to Cell objects '''

        return self.cells

    def get_all_cell_ids(self):
        ''' Returns the id of each cell in the dataset '''

        return self.cells.keys()

    
    def load_from_db(self, db_path):
        '''
        Loads the Dataset instance from a saved file

        db_path (string) : path to database file (.h5)
        '''
        self.db_path = db_path 

        db = h5py.File(db_path, "r")
        cells =  db.get('/cells')
        #Load attributes
        self.voxel_dim = db.attrs['voxel_dim']
        self.volume_dim = db.attrs['volume_dim']
        #Load cells
        for cell_id in cells:
            self.load_cell_from_db(cells, cell_id)

    def load_cell_from_db(self, cells, cell_id):
        '''
        Helper function to load_from_db, loads a cell into the dataset

        cells : database object
        cell_id : the id of the cell to load
        '''

        cell = cells.get(cell_id)
        cell_object = Cell(cell_id,  cell.attrs['cell_type'])

        #Load all regions
        for region_type in cell:
            self.load_region_from_db(cell, cell_object, region_type)
        
        self.cells[str(cell_id)] = cell_object

    def load_region_from_db(self, cell, cell_object, region_type):
        '''
        Helper function to load_cell, loads all annotated cell regions into the dataset

        cell : database object
        cell_object (Cell) : Cell object from which to query the regions
        region_type (string) : the type of cell region annotations to load
        '''

        region_ids = cell.get(region_type)
        for region_id in region_ids:
            #Get data
            region = region_ids.get(region_id)
            voxels = np.array(region['full'], np.uint32)
            membrane =  np.array(region['membrane'], np.uint32)
            #Create Region object and add to cell
            region_object = CellRegion(region_type, region_id, voxels, membrane)
            cell_object.add_region(region_object)

    def load_from_image_stack(self, voxel_dim, cell_stack, region_annotations = [], region_types = [], cell_type_dict = None):
        '''
        Loads a ground truth stack into the dataset object

        db_path (string) : the path to the new db_file
        voxel_dim (tuple) : dimension of a voxel in nm
        cell_stack (numpy Z x X x Y array, uint32) : stack of ground truth data. Each voxel should have a cell_id as value, and 0 if extra-cellular space.
        region_annotations : list of stacks with same shape as cell_stack, where voxels representing certain cell regions are annotated with region ids (for instance a synapse stack)
        region_types : list of strings representing the types of region annotations given. For each stack in region_annotations, region_types stores the corresponding annotated region.
        cell_type_dict (dict {cell_id(int) : type(string)}) : map from cell_id to type (optional)
        '''

        self.voxel_dim = voxel_dim
        self.volume_dim = cell_stack.shape

        #Find all cell_ids in the volume
        non_zero = cell_stack.nonzero()# Returns a (Z, X, Y) tuple
        all_cell_ids = np.unique(cell_stack[non_zero]) #Takes the unique values in the non zero stack

        for cell_id in all_cell_ids:
            #Create Cell object
            cell = Cell(cell_id)
            if cell_type_dict and cell_type_dict.has_key(cell_id): 
                cell.set_cell_type(cell_type_dict[cell_id])
            #Load main voxels and compute membrane
            self.load_cell_from_stack(cell, cell_stack, non_zero)

            #Load other annotations
            for i in range(len(region_annotations)):
                #Load region annotations
                self.load_region_from_stack(cell, region_annotations[i], region_types[i])

            self.add_cell(cell)

    def load_cell_from_stack(self, cell, cell_stack, non_zero):
        '''
        Helper function to load_image_stack, loads a cell by computing its voxels in three dimension as well as its edges (i.e membrane)

        cell (Cell) : the cell object to which we will add the full and membrane locations
        cell_stack (numpy Z x X x Y array, uint32) : stack of ground truth data. Each voxel should have a cell_id as value, and 0 if extra-cellular space.
        non_zero ((Z, X, Y) tuple, uint32) : the non_zero locations, for faster loading 
        '''

        #Find cell voxels
        voxel_indices = np.where(cell_stack[non_zero] == int(cell.get_cell_id()))
        voxel_list = np.transpose(non_zero)[voxel_indices, :][0]
        #tranpose creates an array with 3 axis, 1 st one of size 1, so this squeezes it.

        #Compute membranes

        membrane_indices = np.zeros(voxel_list.shape[0], np.bool)
        for i in range(voxel_list.shape[0]):
            (z, x, y) = voxel_list[i, :]
            membrane_indices[i] = np.unique(cell_stack[z, x - 1 : x + 2, y - 1 : y + 2]).size > 1

        membrane_list = voxel_list[membrane_indices, :]

        #Set the region 'full' of the cell, this method is a shortcut for adding a CellRegion with region_id = 1 and region_type = 'full'
        cell.set_full_cell(voxels = voxel_list, membrane = membrane_list)



    def load_region_from_stack(self, cell, region_stack, region_type):
        '''
        Helper function to load_image_stack. Loads an annotated region. 
        Such region could be synapses, dendrites, axons, mitochondria.. 
        The idea is to use image stacks which have the region_id as value. By superposing the cell_stack,
        we can easily reconstruct which regions belong to what cell and further divide them by region_id, 
        for instance to label a specific synapse and not others.

        cell (Cell) : the cell object to which we will add the region annotation
        region_stack (numpy Z x X x Y array, uint32) : stack of ground truth data. Each voxel should have a region_id as value, and 0 if extra-cellular space.
        region_type (string) : the type of region this data is annotating, for future querying.
        '''

        voxel_list = cell.get_full_cell()
        membrane_list = cell.get_full_cell(membrane_only = True)

        all_region_ids = region_stack[voxel_list.transpose()].unique()[1:] #remove 0
        for region_id in all_region_ids:
            #Find locations
            voxel_indices = np.where(region_stack[voxel_list.transpose()] == region_id)
            membrane_indices = np.where(region_stack[membrane_list.transpose()] == region_id)
            #create Region object and add to cell
            region_voxel_list = voxel_list[voxel_indices, :][0]
            region_membrane_list  = membrane_list[membrane_indices, :][0]
            region = CellRegion(region_type, region_id, region_voxel_list, region_membrane_list)
            cell.add_region(region)


    def save_to_db(self, db_path):
        '''
        Save the Dataset in its current configuration at the given path

        db_path (string) : path to database file (.hdf5)
        '''

        self.db_path = db_path

        db = h5py.File(db_path, "w")

        db.attrs['volume_dim'] = self.volume_dim
        db.attrs['voxel_dim'] = self.voxel_dim
        cells = db.create_group('cells')

        for cell_id in self.cells.keys():
            cell_object = self.cells[cell_id]
            self.save_cell(cells, cell_object)

    def save_cell(self, cells, cell_object):
        '''
        Herlper function to save_to_db, save cell to the given db

        cells : database object
        cell_object (Cell) : Cell object to save
        '''

        cell = cells.create_group(str(cell_object.get_cell_id()))
        cell.attrs['cell_type'] = cell_object.get_cell_type()

        for region_type in cell_object.get_region_types():
            self.save_region(cell, cell_object, region_type)

    def save_region(self, cell, cell_object, region_type):
        '''
        Helper function to save_cell, save all regions of the given cell

        cell : database object for the specific cell
        cell_object (Cell) : Cell object to save
        region_type (string) : the type of cell region annotations to save
        '''

        region_group = cell.create_group(region_type)
        for region_id in cell_object.get_all_region_ids(region_type):
            region = region_group.create_group(str(region_id))
            voxels = cell_object.get_region(region_type, region_id, membrane_only = False)
            membrane = cell_object.get_region(region_type, region_id, membrane_only = True)
            region.create_dataset('full', data = voxels)
            region.create_dataset('membrane', data = membrane)
    
    def get_parameters(self):
        '''
        Returns a dicitonary containing the ground truth parameters
        '''
        params = {}
        params['voxel_dim'] = list(self.voxel_dim) # nm
        params['volume_dim'] = list(self.volume_dim) # number of voxel per dimension
        params['data_path'] = self.db_path # path to the temporary db file

        return params    


class Fluorset: 

    '''
    Interface to access fluorophore data. The data is located in the fluordata folder.
    This is a pre-initialized dataset which is used to query flurophore parameters in the labeling simualtion.

    '''
    def __init__(self):
        '''
        Init method, hard-coded.
        all_fluors (list of FLuorophore objects) : list of all fluorophore instances in fluors.py
        '''

        self.all_fluors = [f.Alexa350(), f.Alexa790(), f.ATTO390(), f.ATTO425(), f.ATTO430LS(), f.ATTO465(), f.ATTO488(), f.ATTO490LS(), f.ATTO550(), f.ATTO647N(), f.ATTO700()]

    def get_all_fluorophores_types(self):
        ''' Returns a list of flurophore types as a list of strings '''

        return [fluor.get_name() for fluor in self.all_fluors]

    def get_fluor(self, name):
        '''
        Returns the desired Fluorophore object. See fluors.py

        name (string) : the name of the fluorophore to query
        '''

        fluors = self.get_all_fluorophores_types()
        index = fluors.index(name)
        return self.all_fluors[index]


