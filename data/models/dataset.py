'''
dataset.py

The Dataset interface for ground truth data and,
the Fluorset interface to access the Fluorophore objects.

'''
import h5py
import os.path
from cell import Cell
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
        #Load attributes
        self.voxel_dim = db.attrs['voxel_dim']
        self.volume_dim = db.attrs['volume_dim']
        #Load cells
        all_cells = db.get('/cells'):
        for cell_id in all_cells:
            cell = db.get('/cells/' + str(cell_id))
            cell_object = Cell(cell_id,  cell.attrs['cell_type'])
            full_region_object = CellRegion(np.array(cell['full'], np.uint32), np.array(cell['membrane'], np.uint32, region_type='full')
            cell_object.add_region(full_region_object)
            #Load regions
            for region_type in db.get('/cells/regions'):
                for region_id in db.get('/cells/regions/' + region_type):
                    region = db.get('/cells/regions/' + region_type + '/' + str(region_id))
                    region_object = CellRegion(np.array(region['full']), np.uint32), np.array(region['membrane'], np.uint32), region_type, region_id)
                    cell_object.add_region(region_object)
            
            self.cells[str(cell_id)] = cell_object

    def save_to_db(self, db_path):
        '''
        Save the Dataset in its current configuration at the given path

        db_path (string) : path to database file (.h5)
        '''

        db = h5py.File(db_path, "w")

        db.attrs['volume_dim'] = self.volume_dim
        db.attrs['voxel_dim'] = self.voxel_dim
        db.create_group('/cells')

        for cell in self.cells:
            cell_id = db.create_group('/cells/' + str(cell.get_cell_id()))
            cell_id.attrs['cell_type'] = cell.get_cell_type()

            data_full = cell.get_full_cell()
            full = self.db.create_dataset('/cells/' + str(cell.get_cell_id()) + '/full', data_full.shape)
            full[...] = data_full

            data_membrane = cell.get_full_cell(membrane_only = True)
            membrane = self.db.create_dataset('/cells/' + str(cell.get_cell_id()) + '/membrane', data_membrane.shape)
            membrane[...] = data_membrane

            for region_type in cell.get_region_types():
                db.create_dataset('/cells/regions/' + region_type)
                for region_id in cell.get_all_region_ids_of_type(region_type):
                    data_full_reg = cell.get_region(region_type, region_id)
                    full_reg = self.db.create_dataset('/cells/regions/' + region_type + '/' + region_id + '/full', data_full_reg)
                    full_reg[...] = data_full_reg

                    data_membrane_reg = cell.get_region(region_type, region_id)
                    membrane_reg = self.db.create_dataset('/cells/regions/' + region_type + '/' + region_id + '/membrane', data_membrane_reg.shape)
                    membrane_reg[...] = data_membrane_reg


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


