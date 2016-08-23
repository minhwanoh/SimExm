'''
dataset.py

The Dataset interface for ground truth data and,
the Fluorset interface to access the Fluorophore objects.

'''
import h5py
import os.path
from cell import Cell, CellRegion
import fluors as f


class Dataset:

    '''
    Interface to ingest, manipulate, and store ground truth datasets.
    Labeling units take Datasets as input.

    Uses hdf5 as storage format.
    Use load_from_db and save_to_db as ways to read from or write to a .h5 file

    '''
   
    def __init__(self, voxel_dim = None, volume_dim = None):
        ''' 
        Init method, sets attributes

        voxel_dim (1x3 int tuple) : the dimensions of a voxel in nm
        volume_dim (1x3 int tuple) : the dimensions of the whole dataset in voxel count
        '''

        self.cells = {}
        self.volume_dim = volume_dim
        self.voxel_dim = voxel_dim
            
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

        db = h5py.File(db_path, "r")
        cells =  db.get('/cells')
        #Load attributes
        self.voxel_dim = db.attrs['voxel_dim']
        self.volume_dim = db.attrs['volume_dim']
        #Load cells
        for cell_id in cells:
            self.load_cell(cells, cell_id)

    def load_cell(self, cells, cell_id):
        '''
        Helper function to load_from_db, loads a cell into the dataset

        cells : database object
        cell_id : the id of the cell to load
        '''

        cell = cells.get(cell_id)
        cell_object = Cell(cell_id,  cell.attrs['cell_type'])

        #Load all regions
        for region_type in cell:
            self.load_region(cell, cell_object, region_type)
        
        self.cells[str(cell_id)] = cell_object

    def load_region(self, cell, cell_object, region_type):
        '''
        Helper function to load_cell, loads all annotated cell regions into the dataset

        cell : database object
        cell_object (Cell) : Cell object from which to query the regions
        region_type (string) : the type of cell region annotations to load
        '''

        region_type = cell.get(region_type)
        for region_id in region_type:
            #Get data
            region = region_type.get(str(region_id))
            voxels = np.array(region['full']), np.uint32)
            membrane =  np.array(region['membrane'], np.uint32)
            #Create Region object and add to cell
            region_object = CellRegion(region_id, voxels, membrane, region_type)
            cell_object.add_region(region_object)


    def save_to_db(self, db_path):
        '''
        Save the Dataset in its current configuration at the given path

        db_path (string) : path to database file (.h5)
        '''

        db = h5py.File(db_path, "w")

        db.attrs['volume_dim'] = self.volume_dim
        db.attrs['voxel_dim'] = self.voxel_dim
        db.create_group('/cells')

        for cell_id in self.cells.keys():
            cell_object = self.cells[cell_id]
            self.save_cell(db, cell_object)

    def save_cell(self, db, cell_object):
        '''
        Herlper function to save_to_db, save cell to the given db

        db : database object
        cell_object (Cell) : Cell object to save
        '''

        cell = db.create_group('/cells/' + str(cell.get_cell_id()))
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

        region_type = cell.create_group(region_type)
        for region_id in cell_object.get_all_region_ids_of_type(region_type):
            region = region_type.create_group(region_id)
            voxels = cell_object.get_region(region_type, region_id, membrane_only = False)
            membrane = cell_object.get_region(region_type, region_id, membrane_only = True)
            region.create_dataset('full', data = voxels)
            region.create_dataset('membrane', data = membrane)
            


class Fluoroset: 

    '''
    Interface to access fluorophore data. The data is located in the fluordata folder.
    This is a pre-initialized dataset which is used to query flurophore parameters in the labeling simualtion.

    '''
    def __init__(self):
        '''
        Init method, hard-coded.
        all_fluors (list of FLuorophore objects) : list of all fluorophore instances in fluors.py
        '''

        self.all_fluors = [f.Alexa350, f.Alexa790, f.ATTO390, f.ATTO425, f.ATTO430LS, f.ATTO465, f.ATTO488, f.ATTO490LS, f.ATTO550, f.ATTO647N, f.ATTO700]

    def get_all_fluorophores_types(self):
        ''' Returns a list of flurophore types as a list of strings '''

        return [fluor.get_name() for fluor in self.all_fluors]

    def get_fluor(self, name):
        '''
        Returns the desired Fluorophore object. See fluors.py

        name (string) : the name of the fluorophore to query
        '''

        all_fluors = get_all_fluorophores_types()
        index = [i.name == name for i in all_fluors].index(1)
        return all_fluors[index]


