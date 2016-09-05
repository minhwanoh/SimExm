'''
sim.py

The Simstack and SimParam class.

'''

from images2gif import writeGif
from tifffile import imsave
from PIL import Image
import numpy as np
import json

class SimStack:
    '''
    Handles simulation outputs. Saves to gif, image_sequence (png) or tiff stack.
    '''

    def __init__(self, image_volume, ground_truth, sim_params):
        '''
        Init method, set attributes.

        image_volume (numpy Z x X x Y x C volume, uint8)  : the  simulation output
        ground_truth (numpy Z x X x Y, uint32 ) : the simulation ground truth output
        sim_params (SimParams) : an instance of the SimParams class, which helps deal with parameter file
        '''

        self.image_sequence = [np.squeeze(image_volume[i, :, :, :]) for i in range(image_volume.shape[0])]
        self.ground_truth = [np.squeeze(ground_truth[i, :, :]) for i in range(ground_truth.shape[0])]
        self.num_channels = image_volume.shape[3]
        self.sim_params = sim_params

        if self.num_channels < 3:
            self.add_empty_channel(num = 3 - self.num_channels)

    def add_empty_channels(self, num = 1):
        '''
        Add an empty channel to each image num times

        num (int) : number of empty channels to add to each image
        '''

        (x, y, c) = self.image_sequence[0].shape
        empty = np.zeros((x, y, num), np.uint8)
        for i in range(len(self.image_sequence)):
            image = self.image_sequence[i]
            new_image = np.concatenate([image, empty], 2)
            self.image_sequence[i] = new_image

    def save_as_gif(self, path, name):
        ''' 
        Save image sequence as gif for immediate visualization 

        path (string) : the destination folder
        name (string) : the name of the file
        '''

        if path[-1] != "/":
            path += "/"
        writeGif(path + name + ".gif", self.image_sequence, duration=0.5)
        writeGif(path + name + "_gt.gif", self.ground_truth, duration=0.5)
        self.sim_params.save(path, name + "_params")

    def save_as_image_sequence(self, path, name):
        '''
        Save image sequence as a sequence of png images in the given directory

        path (string) : the destination folder
        name (string) : the name of the file
        '''

        if path[-1] != "/":
            path += "/"
        for i in range(len(self.image_sequence)):
            im = Image.fromarray(self.image_sequence[i], 'RGB')
            gt = Image.fromarray(self.ground_truth[i], 'L')
            im.save(path + name + "/" + "image_" + str(i) + ".png")
            gt.save(path + name + "_gt/" + "image_" + str(i) + ".png")
        self.sim_params.save(path, name + "_params")

    def save_as_tiff(self, path, name, sixteen_bit_mode = False):
        '''
        Save image sequence as multi-page tiff stack

        path (string) : the destination folder
        name (string) : the name of the file
        sixteen_bit_mode (boolean) : if True, the Tiff stack has 16 bits per channel, default is 8 bits
        '''
        if path[-1] != "/": path += "/"
        if sixteen_bit_mode:
            imsave(path + name + ".tif", np.array(self.image_sequence, np.uint16), photometric='rgb')
        else:
            imsave(path + name + ".tif", np.array(self.image_sequence, np.uint8), photometric='rgb')
        imsave(path + name + "_gt.tif", np.array(self.ground_truth, np.float32), photometric='minisblack')
        self.sim_params.save(path, name + "_params")


class SimParams:

    '''
    A class to store simulation parameters
    '''

    def __init__(self, volume_dim, gt_params, labeling_params, expansion_params, optics_params):
        '''
        Init method, set attributes
            
        volume_dim (tuple (x, y, z ) int): the size of the output sim
        gt_params (dict) : dictionary from parameter to value
        labeling_params (dict) : dictionary from parameter to value
        optics_params (dict) : dictionary from parameter to value
        expansion_params (dict) : dictionary from parameter to value
        '''

        self.volume_dim = volume_dim
        self.gt_params = gt_params
        self.labeling_params = labeling_params
        self.optics_params = optics_params
        self.expansion_params = expansion_params

    def save(self, path, name):

        '''
        Save image sequence as multi-page tiff stack

        path (string) : the destination folder
        name (string) : the name of the file
        '''

        if path[-1] != "/": path += "/"

        pre_expansion_space_voxel_dim = [self.optics_params['pixel_size'] / self.optics_params['objective_factor'],\
          self.optics_params['pixel_size'] / self.optics_params['objective_factor'], self.optics_params['focal_plane_depth']]
        post_expansion_space_voxel_dim = [pre_expansion_space_voxel_dim[i] / self.expansion_params['expansion_factor'] for i in range(3)]

        total_gt_size = [self.gt_params['volume_dim'][i] * self.gt_params['voxel_dim'][i] for i in range(3)]
        total_sim_size_pre = [self.volume_dim[i] * pre_expansion_space_voxel_dim[i] for i in range(3)]
        total_sim_size_post = [self.volume_dim[i] * post_expansion_space_voxel_dim[i] for i in range(3)]

        parameters = { 'volume_dim' : self.volume_dim, \
                    'ground_truth' : self.gt_params, 'labeling' : self.labeling_params, 'optics': self.optics_params, 'expansion' : self.expansion_params, \
                    'pre_expansion_space_voxel_dim' : pre_expansion_space_voxel_dim, 'post_expansion_space_voxel_dim' : post_expansion_space_voxel_dim,
                    'total_ground_truth_size' : total_gt_size, 'total_sim_size_pre_exp_space' : total_sim_size_pre, 'total_sim_size_posT_exp_space' : total_sim_size_post\
                    }
        with open(path + name + '.txt', 'w') as f:
            json.dump(parameters, f)



