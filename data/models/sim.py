'''
sim.py

The Simstack and SimParam class.

'''

from images2gif import writeGif
from tifffile import imsave
import json

class SimStack:
    '''
    Handles simulation outputs. Saves to gif, image_sequence (png) or tiff stack.
    '''

    def __init__(self, image_sequence, ground_truth, num_channels):
        '''
        Init method, set attributes.

        image_sequence : list of numpy X x Y x C, uint8 arrays
        num_channels (int) : the number of channels in each image
        '''

        self.image_sequence = image_sequence
        self.num_channels = num_channels

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

    def save_as_gif(path, name):
        ''' 
        Save image sequence as gif for immediate visualization 

        path (string) : the destination folder
        name (string) : the name of the file
        '''

        if path[-1] != "/":
            path += "/"
        writeGif(path + name + ".gif", self.image_sequence, duration=0.5)
        writeGif(path + name + "_gt.gif", self.ground_truth, duration=0.5)

    def save_as_image_sequence(path, name):
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

    def save_as_tiff_stack(path, name, sixteen_bit_mode = False):
        '''
        Save image sequence as multi-page tiff stack

        path (string) : the destination folder
        name (string) : the name of the file
        sixteen_bit_mode (boolean) : if True, the Tiff stack has 16 bits per channel, default is 8 bits
        '''
        if path[-1] != "/": path += "/"
        if sixteen_bit_mode:
            imsave(path + name + ".tif", np.array(image_sequence, np.uint16))
        else:
            imsave(path + name + ".tif", np.array(image_sequence, np.uint8))
        imsave(path + name + "_gt.tif", np.array(image_sequence, np.float32))


class SimParam:

    '''
    A class to store simulation parameters
    '''

    def __init__(self, labeling_params, optical_params, expansion_params):
        '''
        Init method, set attributes

        labeling_params (dict) : dictionary from parameter to value
        optical_params (dict) : dictionary from parameter to value
        expansion_params (dict) : dictionary from parameter to value
        '''

        self.labeling_params = labeling_params
        self.optical_params = optical_params
        self.expansion_params = expansion_params

    def save_param(self, path, name):

        '''
        Save image sequence as multi-page tiff stack

        path (string) : the destination folder
        name (string) : the name of the file
        '''

        if path[-1] != "/": path += "/"

        labeling_params = labeling_unit.get_param_dict()
        optics_params = optical_unit.get_param_dict()
        parameters = {'labeling_params': self.labeling_params, 'optics_params': self.optics_params, 'expansion_params' = self.expansion_params}
        with open(path + name + '.txt', 'w') as f:
            json.dump(parameters, f)



