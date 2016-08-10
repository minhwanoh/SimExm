# cython: profile=True
import numpy as np
cimport cython
cimport numpy as np
from cython.parallel cimport prange, parallel

############# WRITE #################

""" Add a new data point to the ground truth data. """
cpdef add_data_point(object ds, int dataset_id, int cell_id, int z_offset, int y_offset, int x_1, int x_2, object region_type, int region_id, bint membrane):
    ds.execute('''INSERT INTO connectomics_data (dataset_id, cell_id, z_offset, y_offset, x_1, x_2, region_type_id, region_id, membrane) VALUES\
     (?,?,?,?,?,?, (SELECT type_id FROM region_types WHERE type=?), ?, ?)''', (dataset_id, cell_id, z_offset, y_offset, x_1, x_2, region_type, region_id, membrane,))

""" Add multiple ground truth data points at the same time using executemany """
cpdef add_data_points(object ds, object data):
    ds.executemany('''INSERT INTO connectomics_data (dataset_id, cell_id, z_offset, y_offset, x_1, x_2, region_type_id, region_id, membrane) VALUES\
     (?,?,?,?,?,?, (SELECT type_id FROM region_types WHERE type=?), ?, ?)''', data)

""" Add a new cell to the ground truth data. """
cpdef add_cell(object ds, int dataset_id, int cell_id, object cell_type):
    ds.execute('''INSERT OR IGNORE INTO cells (dataset_id, cell_id, type_id) VALUES (?,?, (SELECT type_id FROM cell_types WHERE type=?))''',(dataset_id, cell_id, cell_type,))

    """ Add a new cell to the ground truth data. """
cpdef add_cells(object ds, object data):
    ds.executemany('''INSERT OR IGNORE INTO cells (dataset_id, cell_id, type_id) VALUES (?,?, (SELECT type_id FROM cell_types WHERE type=?))''', data)


############# READ #################

""" Returns the total number of cells present in the dataset """
cpdef int get_number_of_cells(object ds, int dataset_id):
    ds.execute('''SELECT COUNT(cell_id) FROM cells WHERE dataset_id=?''', (dataset_id,))
    number, = ds.fetchone()
    return <int>number

""" Returns the number of neurons present in the dataset """
cpdef int get_number_of_cells_of_type(object ds, int dataset_id, cell_type):
    ds.execute('''SELECT COUNT(cell_id) FROM cells JOIN cell_types ON cells.type_id=cell_types.type_id \
        WHERE cells.dataset_id=? AND cell_types.type=?''', (dataset_id, cell_type,))
    number, = ds.fetchone()
    return <int>number

""" Returns the ids of all cells present in the dataset """
cpdef np.ndarray[np.int64_t] get_all_cell_ids(object ds, int dataset_id):
    ds.execute('''SELECT cell_id FROM cells WHERE dataset_id=?''', (dataset_id,))
    return np.array([i[0] for i in ds.fetchall()], dtype=(np.int))

""" Returns the ids of all cells of the given type present in the dataset """
cpdef np.ndarray[np.int64_t] get_all_ids_of_type(object ds, int dataset_id, object cell_type):
    ds.execute('''SELECT cell_id FROM cells JOIN cell_types ON cells.type_id=cell_types.type_id \
        WHERE cells.dataset_id=? AND cell_types.type=?''', (dataset_id, cell_type,))
    return np.array([i[0] for i in ds.fetchall()], dtype=(np.int))

""" Returns the cell type of the corresponding cell"""
cpdef get_cell_type(object ds, int dataset_id, int cell_id):
    ds.execute('''SELECT type FROM cell_types JOIN cells ON cells.type_id=cell_types.type_id \
        WHERE cells.dataset_id=? AND cells.cell_id=?''', (dataset_id, cell_id,))
    cell_type, = ds.fetchone()
    return cell_type

""" Returns an array of all voxel locations associated with the given cell in the form of (x, y, z) tuples """
cpdef np.ndarray[np.int64_t, ndim=2] get_all_cell(object ds, int dataset_id, int cell_id):
    ds.execute('''SELECT x_1, x_2, y_offset, z_offset FROM connectomics_data WHERE dataset_id=? AND cell_id=?''', (dataset_id, cell_id,))
    data = ds.fetchall()
    return np.empty(0) if len(data) == 0 else np.concatenate([get_all_points(i[0], i[1], i[2], i[3]) for i in data])

""" Returns the number of voxels, in the volume, which are associated with the given cell id """
cpdef long long get_number_of_voxels_in_cell(object ds, int dataset_id, int cell_id):
    ds.execute('''SELECT SUM(x_2 + 1 - x_1) FROM connectomics_data WHERE dataset_id=? AND cell_id=?''', (dataset_id, cell_id,))
    number_of_voxels, = ds.fetchone()
    return number_of_voxels

""" Returns the number of membrane voxels, in the volume, which are associated with the given cell id """
cpdef long long get_number_of_membrane_voxels_in_cell(object ds, int dataset_id, int cell_id):
    ds.execute('''SELECT SUM(x_2 + 1 - x_1) FROM connectomics_data WHERE dataset_id=? AND cell_id=? AND membrane = 1''', (dataset_id, cell_id,))
    number_of_voxels, = ds.fetchone()
    return number_of_voxels

""" Returns an array of all voxel locations associated with the given cell body in the form of (x, y, z) tuples """
cpdef np.ndarray[np.int64_t, ndim=2] get_cell_body(object ds, int dataset_id, int cell_id):
    ds.execute('''SELECT x_1, x_2, y_offset, z_offset FROM connectomics_data WHERE dataset_id=? AND cell_id=? AND membrane=0''', (dataset_id, cell_id,))
    data = ds.fetchall()
    return np.empty(0) if len(data) == 0 else np.concatenate([get_all_points(i[0], i[1], i[2], i[3]) for i in data])

""" Returns an array of all voxel locations on the membrane of the given cell in the form of (x, y, z) tuples """
cpdef np.ndarray[np.int64_t, ndim=2] get_membrane(object ds, int dataset_id, int cell_id):
    ds.execute('''SELECT x_1, x_2, y_offset, z_offset FROM connectomics_data WHERE dataset_id=? AND cell_id=? AND membrane=1''', (dataset_id, cell_id,))
    data = ds.fetchall()
    return np.empty(0) if len(data) == 0 else np.concatenate([get_all_points(i[0], i[1], i[2], i[3]) for i in data])

""" Returns an array of all voxel locations on the given region of the cell in the form of (x, y, z) tuples """
cpdef np.ndarray[np.int64_t, ndim=2] get_region(object ds, int dataset_id, int cell_id, object region_type):
    ds.execute('''SELECT x_1, x_2, y_offset, z_offset FROM connectomics_data JOIN region_types ON connectomics_data.region_type_id=region_types.type_id WHERE\
     connectomics_data.dataset_id=? AND connectomics_data.cell_id=? AND region_types.type=?''',(dataset_id, cell_id, region_type,))
    data = ds.fetchall()
    return np.empty(0) if len(data) == 0 else np.concatenate([get_all_points(i[0], i[1], i[2], i[3]) for i in data])

""" Returns an array of all voxel locations on the given region of the cell in the form of (x, y, z) tuples """
cpdef np.ndarray[np.int64_t, ndim=2] get_local_region(object ds, int dataset_id, int cell_id, object region_type, int region_id):
    ds.execute('''SELECT x_1, x_2, y_offset, z_offset FROM connectomics_data JOIN region_types ON connectomics_data.region_type_id=region_types.type_id WHERE\
     connectomics_data.dataset_id=? AND connectomics_data.cell_id=? AND region_types.type=? AND connectomics_data.region_id=region_id''',(dataset_id, cell_id, region_type, region_id))
    data = ds.fetchall()
    return np.empty(0) if len(data) == 0 else np.concatenate([get_all_points(i[0], i[1], i[2], i[3]) for i in data])

############# REAL_SPACE #################

""" Returns the cell_id of the cell that the voxel at location (x, y, z) belongs to. 0 indicates extracellular space """
cpdef int get_cell_id(object ds, int dataset_id, int x, int y, int z):
    ds.execute('''SELECT cell_id FROM connectomics_data WHERE dataset_id=? AND z=? AND y=? AND ? BETWEEN x_1 AND x_2''',(dataset_id, z, y, x))
    cell_id, = ds.fetchone()
    return cell_id if cell_id is not None else 0

cpdef np.ndarray[np.uint32_t] get_cell_ids_in_volume(object ds, int dataset_id, object bounds):
    (x1, x2, y1, y2, z1, z2) = bounds
    ds.execute('''SELECT UNIQUE cell_id FROM connectomics_data WHERE dataset_id=?\
     AND (z_offset BETWEEN ? AND ?) AND (y_offset BETWEEN ? AND ?) AND ((x_1 BETWEEN  ? AND ?) OR (x_2 BETWEEN ? AND ?))''',\
      (y1, y2 - 1, dataset_id, z1, z2 - 1, x1, x2 - 1, y1, y2 - 1, y1, y2 - 1,))
    return np.array([i[0] for i in ds.fetchall()], dtype=(np.uint32))


""" Returns a dictionary with keys: cell_id and values:[array of membrane voxels], restricted by the given bounds"""
cpdef object get_voxels_in_volume_by_cell(object ds, int dataset_id, object bounds, bint include_body, bint include_membrane):
    cdef int x1, x2, y1, y2, z2, i, length
    cdef np.ndarray[np.int32_t, ndim=2] data, cell_block
    (x1, x2, y1, y2, z1, z2) = bounds
    cell_dict = {}

    if include_membrane and include_body:
        append_string = ''''''
    elif include_membrane and not include_body:
        append_string = ''' AND membrane = 1 '''
    elif not include_membrane and include_body:
        append_string = ''' AND membrane = 0 '''
    else: 
        return cell_dict

    ds.execute('''SELECT cell_id, max(x_1, ?), min(x_2, ?), y_offset, z_offset FROM connectomics_data WHERE dataset_id=?\
     AND (z_offset BETWEEN ? AND ?) AND (y_offset BETWEEN ? AND ?) AND ((x_1 BETWEEN  ? AND ?) OR (x_2 BETWEEN ? AND ?))''' + append_string + ''' ORDER BY cell_id''',\
      (x1, x2 - 1, dataset_id, z1, z2 - 1, y1, y2 - 1, x1, x2 - 1, x1, x2 - 1,))

    datab = ds.fetchall()
    if len(datab) == 0:
        return cell_dict
    else:
        data = np.array(datab, np.int32)
    splitted_data = np.split(data, np.diff(data[:,0]).nonzero()[0] + 1)
    lenght = len(splitted_data)
    for i in range(lenght):
        cell_block = splitted_data[i]
        cell_id = cell_block[0][0]
        cell_dict[cell_id] = np.concatenate([get_all_points(cell_data[1] - x1, cell_data[2] - x1, cell_data[3] - y1, cell_data[4] - z1) for cell_data in cell_block])
 
    return cell_dict


""" Returns a dictionary with keys: cell_id and values:[array of membrane voxels], restricted by the given bounds"""
cpdef np.ndarray[np.uint32_t, ndim=3] get_membranes_in_volume(object ds, int dataset_id, object bounds):
    cdef int x1, x2, y1, y2, z2, x, y, z, length, cell_id
    cdef np.ndarray[np.int16_t, ndim=2] locations
    (x1, x2, y1, y2, z1, z2) = bounds
    cdef np.ndarray[np.uint32_t, ndim=3] chunk = np.zeros((x2 - x1, y2 - y1, z2 - z1), np.uint32)
    ds.execute('''SELECT cell_id, max(x_1, ?), min(x_2, ?), y_offset, z_offset FROM connectomics_data WHERE dataset_id=?\
     AND (z_offset BETWEEN ? AND ?) AND (y_offset BETWEEN ? AND ?) AND ((x_1 BETWEEN  ? AND ?) OR (x_2 BETWEEN ? AND ?)) AND membrane = 1''',\
      (x1, x2 - 1, dataset_id, z1, z2 - 1, y1, y2 - 1, x1, x2 - 1, x1, x2 - 1,))
    data = ds.fetchall()
    for data_point in data:
        cell_id = data_point[0]
        locations = get_all_points(data_point[1] - x1, data_point[2] - x1, data_point[3] - y1, data_point[4] - z1)
        for (x, y, z) in locations:
            chunk[x, y, z] = cell_id
    return chunk

""" Returns a dictionary with keys: cell_id and values:[array of membrane voxels], restricted by the given bounds"""
cpdef np.ndarray[np.uint32_t, ndim=3] get_bodies_in_volume(object ds, int dataset_id, object bounds):
    cdef int x1, x2, y1, y2, z2, x, y, z, length, cell_id
    cdef np.ndarray[np.int16_t, ndim=2] locations
    (x1, x2, y1, y2, z1, z2) = bounds
    cdef np.ndarray[np.uint32_t, ndim=3] chunk = np.zeros((x2 - x1, y2 - y1, z2 - z1), np.uint32)
    ds.execute('''SELECT cell_id, max(x_1, ?), min(x_2, ?), y_offset, z_offset FROM connectomics_data WHERE dataset_id=?\
     AND (z_offset BETWEEN ? AND ?) AND (y_offset BETWEEN ? AND ?) AND ((x_1 BETWEEN  ? AND ?) OR (x_2 BETWEEN ? AND ?))''',\
      (x1, x2 - 1, dataset_id, z1, z2 - 1, y1, y2 - 1, x1, x2 - 1, x1, x2 - 1,))
    data = ds.fetchall()
    for data_point in data:
        cell_id = data_point[0]
        locations = get_all_points(data_point[1] - x1, data_point[2] - x1, data_point[3] - y1, data_point[4] - z1)
        for (x, y, z) in locations:
            chunk[x, y, z] = cell_id
    return chunk

""" Returns the membrane_weighting for the given dataset. This is a number used in the simulation. """
cpdef float get_membrane_weighting(object ds, int dataset_id):
    ds.execute('''SELECT SUM(x_2 + 1 - x_1) FROM connectomics_data WHERE dataset_id = ? AND membrane = 1 GROUP BY cell_id''', (dataset_id,))
    membrane = ds.fetchall()
    ds.execute('''SELECT SUM(x_2 + 1 - x_1) FROM connectomics_data WHERE dataset_id = ? GROUP BY cell_id''', (dataset_id,))
    all_voxels = ds.fetchall()
    ratios = np.array(membrane[::], np.float) / np.array(all_voxels[::], np.float) ** (2.0 / 3.0)
    return np.mean(ratios)

############## HELPERS ###############

""" A helper function to unpack the compressed data """
@cython.boundscheck(False)
cdef long[:,:] get_all_points(int x_1, int x_2, int y_offset, int z_offset):
    cdef long[:,:] buffer = np.zeros((x_2 - x_1 + 1, 3), np.int64)
    cdef int index, loc
    with nogil, parallel(num_threads=1):
        for loc in prange(x_1, x_2 + 1, schedule='dynamic'):
            index = loc - x_1
            buffer[index,0] = loc
            buffer[index,1] = y_offset
            buffer[index,2] = z_offset
    return buffer




