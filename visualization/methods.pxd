import cython
import numpy as np
cimport numpy as np

### Visualization ###

''' Create composite from list of neurons in a given volume '''
cpdef object get_composite(object ds, object cell_dict, object dim, max_number_cell = ?)

''' Save image sequence as gif for immediate visualization '''
cpdef make_gif(object image_sequence, object path, object name)

''' Save image sequence as a sequence of png images in the given directory '''
cpdef save_image_sequence(object image_sequence, object path, object name, rgb)

''' Save image sequence as multi-page tiff stack '''
cpdef make_tiff_stack(object image_sequence, object path, object name, rgb, sixteen_bit_mode)
