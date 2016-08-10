import sqlite3
import os

class DatabaseSession:
    def __init__(self, conn_string):
        initalized = os.path.isfile(conn_string)
        self.path = os.path.abspath(os.path.dirname(__file__))
        self.database_connection = sqlite3.connect(conn_string)
        self.database_connection.isolation_level = None
        self.cursor = self.database_connection.cursor()
        if not initalized:
            self.build_database()
            self.build_indexing()
        
    def begin(self):
        self.cursor.execute('''BEGIN''')

    def commit(self):
        self.cursor.execute('''COMMIT''')

    def get_path(self):
        return self.path

    def execute(self, query, values=None):
        if values is not None:
            return self.cursor.execute(query, values)
        return self.cursor.execute(query)

    def executemany(self, query, values=None):
        if values is not None:
            return self.cursor.executemany(query, values)
        return self.cursor.executemany(query)

    def fetchall(self):
        return self.cursor.fetchall()

    def fetchone(self):
        return self.cursor.fetchone()

    def close(self):
        self.database_connection.close()

    def build_database(self):
        self.begin()

        self.execute('''CREATE TABLE datasets (dataset_id INTEGER PRIMARY KEY AUTOINCREMENT, name TEXT, type_id INTEGER)''')

        self.execute('''CREATE TABLE ground_truth_datasets (dataset_id INTEGER PRIMARY KEY, name TEXT,\
         volume_width INTEGER, volume_height INTEGER, volume_depth INTEGER, voxel_width_nm INTEGER,\
          voxel_height_nm INTEGER, voxel_depth_nm INTEGER, source TEXT, reference TEXT, comments TEXT)''')

        self.execute('''CREATE TABLE exm_datasets (dataset_id INTEGER PRIMARY KEY, name TEXT, \
         expansion_factor INTEGER, labeling_type_id INTEGER, optical_type_id INTEGER, volume_width INTEGER, volume_height INTEGER, volume_depth INTEGER, voxel_width_nm INTEGER,\
          voxel_height_nm INTEGER, voxel_depth_nm INTEGER, source TEXT, reference TEXT, comments TEXT)''')

        self.execute('''CREATE TABLE simulation_datasets (dataset_id INTEGER PRIMARY KEY, name TEXT, original_dataset_id INTEGER, \
         expansion_factor INTEGER, labeling_type_id INTEGER, optical_type_id INTEGER, parameters TEXT, volume_width INTEGER, volume_height INTEGER, volume_depth INTEGER, voxel_width_nm INTEGER,\
          voxel_height_nm INTEGER, voxel_depth_nm INTEGER)''')

        self.execute('''CREATE TABLE segmentation_datasets (dataset_id INTEGER PRIMARY KEY, original_dataset_id INTEGER,\
         segmentation_method INTEGER, model_id INTEGER, parameters TEXT)''')

        #Connectomics
        self.execute('''CREATE TABLE connectomics_data (dataset_id INTEGER, cell_id INTEGER,\
         z_offset INTEGER, y_offset INTEGER, x_1 INTEGER, x_2 INTEGER, region_type_id INTEGER, region_id INTEGER, membrane BOOLEAN)''')
        self.execute('''CREATE TABLE cells (dataset_id INTEGER, cell_id INTEGER, type_id INTEGER, comments TEXT)''')

        #Exm
        self.execute('''CREATE TABLE exm_data (dataset_id INTEGER, z_offset INTEGER, channel INTEGER, image_path TEXT)''')

        ### Neural net models
        self.execute('''CREATE TABLE neural_net_models (model_id INTEGER PRIMARY KEY, training_dataset_id INTEGER, label_dataset_id INTEGER, batch_size INTEGER, validation_size INTEGER, load_path TEXT)''')


        ### Enumerations
        self.execute('''CREATE TABLE dataset_types (type_id INTEGER PRIMARY KEY AUTOINCREMENT, type TEXT)''')
        self.execute('''CREATE TABLE region_types (type_id INTEGER PRIMARY KEY AUTOINCREMENT, type TEXT)''')
        self.execute('''CREATE TABLE cell_types (type_id INTEGER PRIMARY KEY AUTOINCREMENT, type TEXT)''')
        self.execute('''CREATE TABLE labeling_methods (type_id INTEGER PRIMARY KEY AUTOINCREMENT, type TEXT)''')
        self.execute('''CREATE TABLE optics_methods (type_id INTEGER PRIMARY KEY AUTOINCREMENT, type TEXT)''')
        self.execute('''CREATE TABLE segmentation_methods (type_id INTEGER PRIMARY KEY AUTOINCREMENT, type TEXT)''')

        dataset_types = ['ground_truth', 'simulation', 'exm', 'segmentation']
        cell_types = ["unknown", "neuron", "glia"]
        region_types = ["unknown", "soma", "axon", "dendrite", "synapse", "mitochondria"]
        labeling_methods = ["Brainbow"]
        optics_methods = ["Confocal"]
        segmentation_methods = ["NeuralNet"]

        for i in xrange(0, len(dataset_types)):
            self.execute('''INSERT INTO dataset_types (type) VALUES (?)''', (dataset_types[i],))

        for i in xrange(0, len(cell_types)):
            self.execute('''INSERT INTO cell_types (type) VALUES (?)''', (cell_types[i],))

        for i in xrange(0, len(region_types)):
            self.execute('''INSERT INTO region_types (type) VALUES (?)''', (region_types[i],))

        for i in xrange(0, len(labeling_methods)):
            self.execute('''INSERT INTO labeling_methods (type) VALUES (?)''', (labeling_methods[i],))

        for i in xrange(0, len(optics_methods)):
            self.execute('''INSERT INTO optics_methods (type) VALUES (?)''', (optics_methods[i],))

        for i in xrange(0, len(segmentation_methods)):
            self.execute('''INSERT INTO segmentation_methods (type) VALUES (?)''', (segmentation_methods[i],))

        self.commit()

    def build_indexing(self):
        self.begin()

        #Ground truth
        self.execute('''CREATE INDEX cell on connectomics_data(dataset_id, cell_id)''')
        self.execute('''CREATE INDEX real_space on connectomics_data(dataset_id, z_offset, y_offset, x_1, x_2)''')

        #Images
        self.execute('''CREATE INDEX slice on exm_data(dataset_id, z_offset, channel)''')

        self.commit()

    def __del__(self):
        if self.database_connection is not None:
            self.database_connection.close()


