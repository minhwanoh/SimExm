import numpy as np
import cython
cimport numpy as np
from PIL import Image
from random import randint
from images2gif import writeGif
from tifffile import imsave

"""
Set of methods to save simulation outputs

"""

''' Create composite from list of neurons in a given volume '''
cpdef object get_composite(object ds, object cell_dict, object dim, max_number_cell = None):
    cell_ids = cell_dict.keys()
    (x, y, z) = dim
    volume = np.zeros((z, x, y, 3), np.uint8)
    bound = min(len(cell_ids), max_number_cell) if max_number_cell is not None else len(cell_ids)
    for i in range(bound):
        cell_id = cell_ids[i]
        (Xsm, Ysm, Zsm) = np.split(cell_dict[cell_id], 3, axis=1)
        colorm = [randint(0, 255) for i in xrange(3)]
        volume[Zsm, Xsm, Ysm, :] = colorm
    return volume

''' Save image sequence as gif for immediate visualization '''
cpdef make_gif(object image_sequence, object path, object name):
    if path[-1] != "/": path += "/"
    writeGif(path + name + ".gif", image_sequence, duration=0.5)

''' Save image sequence as a sequence of png images in the given directory '''
cpdef save_image_sequence(object image_sequence, object path, object name, rgb):
    if path[-1] != "/": path += "/"
    if rgb:
        imtype = 'RGB'
    else:
        imtype = 'L'
    for i in xrange(len(image_sequence)):
        im = Image.fromarray(image_sequence[i], imtype)
        im.save(path + name + "/" + "image_" + str(i) + ".png")

''' Save image sequence as multi-page tiff stack '''
cpdef make_tiff_stack(object image_sequence, object path, object name, rgb, sixteen_bit_mode):
    if path[-1] != "/": path += "/"
    if rgb:
        if sixteen_bit_mode:
            imsave(path + name + ".tif", np.array(image_sequence, np.uint16))
        else:
            imsave(path + name + ".tif", np.array(image_sequence, np.uint8))
    else:
        imsave(path + name + ".tif", np.array(image_sequence, np.float32))
