import cython
import numpy as np
cimport numpy as np

### GROUND_TRUTH

@cython.boundscheck(False)
cpdef object ingest_slice(object ds, int dataset_id, int z_offset, np.ndarray[np.uint32_t, ndim=2] segment_slice, np.ndarray[np.uint16_t, ndim=2] synapse_slice,\
 np.ndarray[np.uint16_t, ndim=2] mito_slice, object cell_type_func, object region_type_func, object type_map):
    cdef np.ndarray[np.uint64_t, ndim=1] data_buffer
    cdef np.ndarray[np.int32_t, ndim=2] non_zero = np.array(segment_slice.nonzero(), np.int32)
    cdef int buffer_size = non_zero[0].size
    data_buffer = np.zeros(buffer_size, np.uint64)
    fill_buffer(data_buffer, segment_slice, synapse_slice, mito_slice, non_zero, buffer_size)
    return generate_data_points(non_zero, data_buffer, dataset_id, buffer_size, z_offset, cell_type_func, region_type_func, type_map)

@cython.boundscheck(False)
cdef object generate_data_points(int[:,:] non_zero, unsigned long long[:] data_buffer, int dataset_id, int buffer_size, int z_offset, object cell_type_func, object region_type_func, object type_map):
   #Get data points
    cdef int i, loc_size
    cdef unsigned long long value, cell_id, membrane, region_id, y_offset, is_synapse
    cdef unsigned long long short_mask = 0x0000000000003FFF #for 14 bit number
    cdef unsigned long long long_mask = 0x00000000FFFFFFFF# For 32 bit number
    cdef np.ndarray[np.uint64_t, ndim=1] loc_1, loc_2
    data_points = []
    cell_points = []
    cell_ids_visited = {}

    loc_2 = np.union1d( np.where(np.diff(data_buffer) != 0)[0], np.where(np.diff(non_zero[1]) != 1)[0] ).astype(np.uint64)
    loc_1 = loc_2[::] + 1
    if buffer_size > 0:
        loc_1 = np.insert(loc_1, 0, 0)
        loc_2 = np.append(loc_2, np.uint64(buffer_size - 1))
    loc_size = loc_1.size

    for i in range(loc_size):
        value = data_buffer[loc_1[i]]
        (y_offset, cell_id, region_id, membrane, is_synapse) = (value >> 48, (value >> 16) & long_mask, (value >> 2) & short_mask, (value >> 1) & 0x1, value & 0x1)
        region_type = region_type_func(cell_id, region_id, is_synapse, type_map)
        cell_type = cell_type_func(cell_id, type_map)
        data_points.append((dataset_id, cell_id, z_offset, y_offset, non_zero[1, loc_1[i]], non_zero[1, loc_2[i]], region_type, region_id, membrane))
        if not cell_ids_visited.has_key(cell_id):
            cell_ids_visited[cell_id] = 1
            cell_points.append((dataset_id, cell_id, cell_type))

    return (data_points, cell_points)

@cython.boundscheck(False)
cdef inline bint is_membrane(unsigned int[:,:] segment_slice, int index_i, int index_j) nogil:
    cdef int u, v
    cdef unsigned int current_cell_id = segment_slice[index_i, index_j]
    for u in range(index_i - 1, index_i + 2):
        for v in range(index_j - 1, index_j + 2):
            if (segment_slice[u, v] != current_cell_id):
                return True
    return False

@cython.boundscheck(False)
cdef inline void fill_buffer(unsigned long long[:] data_buffer, unsigned int[:,:] segment_slice, unsigned short[:,:] synapse_slice,\
 unsigned short[:,:] mito_slice, int[:,:] non_zero, int buffer_size) nogil:
    cdef int i, index_i, index_j
    cdef unsigned short region_id, mito_id, synapse_id
    for i in range(buffer_size):
        index_i = non_zero[0, i]
        index_j = non_zero[1, i]
        mito_id = mito_slice[index_i, index_j]
        synapse_id = synapse_slice[index_i, index_j]
        region_id = mito_id if synapse_id == 0 else synapse_id
        data_buffer[i] = ((<unsigned long long > index_i) << 48) |\
         ((<unsigned long long> segment_slice[index_i, index_j]) << 16) |\
         ((<unsigned long long> region_id) << 2) |\
         ((<unsigned long long> is_membrane(segment_slice, index_i, index_j)) << 1) |\
         (<unsigned long long> (synapse_id != 0))