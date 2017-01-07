# Provided under the MIT License (MIT)
# Copyright (c) 2016 Jeremy Wohlwend

# Permission is hereby granted, free of charge, to any person obtaining 
# a copy of this software and associated documentation files (the "Software"),
# to deal in the Software without restriction, including without limitation
# the rights to use, copy, modify, merge, publish, distribute, sublicense,
# and/or sell copies of the Software, and to permit persons to whom the Software
# is furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included
# in all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
# INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
# IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
# DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
# TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
# OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

"""
run.py

Main script to run the simulation.
"""

import sys
from configobj import ConfigObj, flatten_errors
from validate import Validator
from src.load import load_gt 
from src.labeling import label
from src.optics import resolve
from src.output import save, save_gt

def run(config):
    """
    Runs the simulation using the given config. 
    For more information on the available parameters, see configspecs.ini and the readme 

    Args:
        config: dicitonary
            the configuration dict outputed by the configobj library
    """
    gt_params = config['groundtruth']
    volume_dim = gt_params['bounds']
    voxel_dim = gt_params['voxel_dim']
    #Remove voxel dim so that we can pass gt_params to the load_gt function
    del gt_params['voxel_dim']

    print "Loading data..."
    gt_dataset = load_gt(**gt_params)

    print "Labeling..."
    labeling_params = config['labeling']
    labeled_volumes, labeled_cells = label(gt_dataset, volume_dim, voxel_dim, labeling_params)

    print "Imaging..."
    expansion_params = config['expansion']
    optics_params = config['optics']
    volumes = resolve(labeled_volumes, volume_dim, voxel_dim, expansion_params, optics_params)
    print "Saving..."
    #Save to desired output
    output_params = config['output']
    save(volumes, **output_params)
    save_gt(gt_dataset, labeled_cells, volume_dim, voxel_dim, expansion_params, optics_params, **output_params)
    print "Done!"

if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise IOError, 'Please give a single config file as argument'
    #Read config file
    config_file = sys.argv[1]
    config = ConfigObj(config_file, list_values = True,  configspec='configspecs.ini')
    #Validate input configuration
    validator = Validator()
    results = config.validate(validator)

    if results != True:
        for (section_list, key, _) in flatten_errors(config, results):
            if key is not None:
                print 'The "%s" key in the section "%s" failed validation' % (key, ', '.join(section_list))
            else:
                print 'The following section was missing:%s ' % ', '.join(section_list)
    #Run simulation
    run(config)