import numpy as np
import os
############# CREATE #################

""" Creates a new ground truth dataset. The new dataset_id is automatically generated. """
def create_new_dataset(ds, name, type):
    ds.execute('''INSERT INTO datasets (name, type_id) VALUES (?, (SELECT type_id FROM dataset_types WHERE type=?))''', (name, type,))
    ds.execute('''SELECT COUNT(dataset_id) FROM datasets''')
    dataset_id, = ds.fetchone()
    return int(dataset_id)

""" Creates a new ground truth dataset. The new dataset_id is automatically generated. """
def create_ground_truth_dataset(ds, name, volume_dimension, voxel_dimension_nm, source, reference, comments):
    (x, y, z) = volume_dimension
    (xnm, ynm, znm) = voxel_dimension_nm
    dataset_id = create_new_dataset(ds, name, 'ground_truth')
    ds.execute('''INSERT INTO ground_truth_datasets (dataset_id, name, volume_width, volume_height, volume_depth, voxel_width_nm, voxel_height_nm,\
        voxel_depth_nm, source, reference, comments) VALUES (?,?,?,?,?,?,?,?,?,?,?)''', (dataset_id, name, x, y, z, xnm, ynm, znm, source, reference, comments,))
    return dataset_id

""" Creates a new ground truth dataset. The new dataset_id is automatically generated. """
def create_simulation_dataset(ds, name, expansion_factor, labeling_method, optical_method, parameters, original_dataset_id, bounds, comments):
    (x1, x2, y1, y2, z1, z2) = bounds
    dataset_id = create_new_dataset(ds, name, 'simulation')
    os.makedirs(ds.get_path() + '/simulation/dataset_' + str(dataset_id))
    ds.execute('''INSERT INTO simulation_datasets (dataset_id, name, expansion_factor, labeling_type_id, optical_type_id,\
      parameters, original_dataset_id, x1, x2, y1, y2, z1, z2, comments) VALUES (?,?,?,(SELECT type_id FROM labeling_methods WHERE type = ?),\
      (SELECT type_id FROM optical_methods WHERE type = ?),?,?,?,?,?,?,?,?,?)''', (dataset_id, name, expansion_factor,\
       labeling_type_id, optical_type_id, parameters, original_dataset_id, x1, x2, y1, y2, z1, z2, comments,))
    return dataset_id

""" Creates a new ground truth dataset. The new dataset_id is automatically generated. """
def create_exm_dataset(ds, name, expansion_factor, labeling_method, optical_method, volume_dimension, voxel_dimension_nm, source, reference, comments):
    (x, y, z) = volume_dimension
    (xnm, ynm, znm) = voxel_dimension_nm
    dataset_id = create_new_dataset(ds, name, 'exm')
    ds.execute('''INSERT INTO exm_datasets (dataset_id, name, expansion_factor, labeling_method_id, optical_method_id, volume_width, volume_height, volume_depth, voxel_width_nm, voxel_height_nm,\
        voxel_depth_nm, source, reference, comments) VALUES (?,?,?,(SELECT type_id FROM labeling_methods WHERE type = ?),\
      (SELECT type_id FROM optical_methods WHERE type = ?),?,?,?,?,?,?,?,?,?)''', (dataset_id, name, expansion_factor, labeling_method, optical_method, x, y, z, xnm, ynm, znm, source, reference, comments,))
    return dataset_id

""" Creates a new ground truth dataset. The new dataset_id is automatically generated. """
def create_segmentation_dataset(ds, name, volume_dimension, voxel_dimension_nm, source, reference, comments):
    (x, y, z) = volume_dimension
    (xnm, ynm, znm) = voxel_dimension_nm
    dataset_id = create_new_dataset(ds, name, 'segmentation')
    ds.execute('''INSERT INTO segmentation_datasets (name, volume_width, volume_height, volume_depth, voxel_width_nm, voxel_height_nm,\
        voxel_depth_nm, source, reference, comments) VALUES (?,?,?,?,?,?,?,?,?,?)''', (name, x, y, z, xnm, ynm, znm, source, reference, comments,))
    return dataset_id

############# INFO #################

""" Returns an array containing the name of all the datasets. """
def get_all_dataset_names_of_type(ds, type):
    ds.execute('''SELECT name FROM datasets WHERE type_id = (SELECT type_id FROM dataset_types WHERE type = ?)''', (type, ))
    return [i[0] for i in ds.fetchall()]
    
""" Returns an array containing the id of all the datasets. """
def get_all_dataset_ids_of_type(ds, type):
    ds.execute('''SELECT dataset_id FROM datasets WHERE type_id = (SELECT type_id FROM dataset_types WHERE type = ?)''', (type, ))
    return np.array([i[0] for i in ds.fetchall()], np.int)

""" Returns the name of the dataset. """
def get_name(ds, dataset_id):
    ds.execute('''SELECT name FROM datasets WHERE dataset_id=?''', (dataset_id,))
    name, = ds.fetchone()
    return name

""" Returns the name of the dataset. """
def get_type(ds, dataset_id):
    ds.execute('''SELECT type from dataset_types WHERE type_id = (SELECT type_id FROM datasets WHERE dataset_id=?)''', (dataset_id,))
    type, = ds.fetchone()
    return type

""" Returns the dataset_id given the name of the dataset """
def  get_dataset_id(ds, name):
    ds.execute('''SELECT dataset_id FROM datasets WHERE name=?''', (name,))
    dataset_id, = ds.fetchone()
    return int(dataset_id)

""" Returns a tuple (x, y, z) representing the size of the dataset in voxels """
def get_volume_dimensions(ds, dataset_id):
    type = get_type(ds, dataset_id)
    query = "SELECT volume_width, volume_height, volume_depth FROM " + type + "_datasets WHERE dataset_id=?"
    ds.execute(query, (dataset_id,))
    (x, y, z) = ds.fetchall()[0]
    return (x, y, z)

""" Returns a tuple (x, y, z) representing the size of a voxel in nm"""
def get_voxel_dimensions_nm(ds, dataset_id):
    type = get_type(ds, dataset_id)
    query = "SELECT voxel_width_nm, voxel_height_nm, voxel_depth_nm FROM " + type + "_datasets WHERE dataset_id=?"
    ds.execute(query, (dataset_id,))
    (x, y, z) = ds.fetchall()[0]
    return (x, y, z)

""" Returns a tuple (x, y, z) representing the size of the dataset in nm """
def get_volume_dimensions_nm(ds, dataset_id):
    v = get_volume_dimensions(ds, dataset_id)
    vn = get_voxel_dimensions_nm(ds, dataset_id)
    x = v[0] * vn[0]
    y = v[1] * vn[1]
    z = v[2] * vn[2]
    return (x, y, z)