import cython
import numpy as np
cimport numpy as np

### CONNECTOMICS ###

cpdef add_data_point(object ds, int dataset_id, int cell_id, int z_offset, int y_offset, int x_1, int x_2, object region_type, int region_id, bint membrane)

cpdef add_data_points(object ds, object data)

cpdef add_cell(object ds, int dataset_id, int cell_id, object cell_type)

cpdef add_cells(object ds, object data)

cpdef int get_number_of_cells(object ds, int dataset_id)

cpdef int get_number_of_cells_of_type(object ds, int dataset_id, object cell_type)

cpdef np.ndarray[np.int64_t] get_all_cell_ids(object ds, int dataset_id)

cpdef np.ndarray[np.int64_t] get_all_ids_of_type(object ds, int dataset_id, object cell_type)

cpdef get_cell_type(object ds, int dataset_id, int cell_id)

cpdef np.ndarray[np.int64_t, ndim=2] get_all_cell(object ds, int dataset_id, int cell_id)

cpdef long long get_number_of_voxels_in_cell(object ds, int dataset_id, int cell_id)

cpdef long long get_number_of_membrane_voxels_in_cell(object ds, int dataset_id, int cell_id)

cpdef np.ndarray[np.int64_t, ndim=2] get_cell_body(object ds, int dataset_id, int cell_id)

cpdef np.ndarray[np.int64_t, ndim=2] get_all_cell(object ds, int dataset_id, int cell_id)

cpdef np.ndarray[np.int64_t, ndim=2] get_region(object ds, int dataset_id, int cell_id, region_type)

cpdef np.ndarray[np.int64_t, ndim=2] get_membrane(object ds, int dataset_id, int cell_id)

cpdef object get_voxels_in_volume_by_cell(object ds, int dataset_id, object bounds, bint include_body, bint include_membrane)

cpdef np.ndarray[np.uint32_t, ndim=3] get_membranes_in_volume(object ds, int dataset_id, object bounds)

cpdef np.ndarray[np.uint32_t, ndim=3] get_bodies_in_volume(object ds, int dataset_id, object bounds)

cpdef float get_membrane_weighting(object ds, int dataset_id)

cpdef int get_cell_id(object ds, int dataset_id, int x, int y, int z)

cdef long[:,:] get_all_points(int x_1, int x_2, int y_offset, int z_offset)
