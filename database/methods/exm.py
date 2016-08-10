import numpy as np
import tifffile as t

############# CREATE #################

""" Add a new data point to the ground truth data. """
def add_image(ds, dataset_id, z_offset, channel, image_data):
    image_path = ds.get_path() + '/simulation/dataset_' + str(dataset_id) +'/' + 'image_' + str(z) + '_channel_' + str(channel) + '.tiff'
    t.imsave(image_path, image_data, compress = 5)
    ds.execute('''INSERT INTO exm_data (dataset_id, z_offset, channel, image_path) VALUES\
     (?,?,?,?)''', (dataset_id, z_offset, channel, image_path,))

""" Add multiple ground truth data points at the same time using executemany """
def add_image_from_path(ds, dataset_id, z_offset,  image_path, channel = None):
    if channel is not None:
        ds.execute('''INSERT INTO exm_data (dataset_id, z_offset, channel, image_path) VALUES\
            (?,?,?,?)''', (dataset_id, z_offset, channel, image_path,))
    else:
        channels = get_number_of_channels(image_path)
        for channel in xrange(len(channels)):
            ds.execute('''INSERT INTO exm_data (dataset_id, z_offset, channel, image_path) VALUES\
            (?,?,?,?)''', (dataset_id, z_offset, channel, image_path,))

def get_number_of_channels(image_path):
    im = t.imread(image_path)
    if len(im.shape) == 2:
        return 1
    else:
        return min(im.shape)






