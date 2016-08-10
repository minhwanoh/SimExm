cimport numpy as np

cpdef object ingest_slice(object ds, int dataset_id, int z_offset, np.ndarray[np.uint32_t, ndim=2] segment_slice, np.ndarray[np.uint16_t, ndim=2] synapse_slice,\
 np.ndarray[np.uint16_t, ndim=2] mito_slice, object cell_type_func, object region_type_func, object true_map)