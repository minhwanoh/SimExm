# Provided under BSD license
# Copyright (c) 2017, Jeremy Wohlwend
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without modification,
# are permitted provided that the following conditions are met:
#
#  - Redistributions of source code must retain the above copyright notice,
#    this list of conditions and the following disclaimer.
#  - Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#  - Neither the name of SimExm nor the names of its contributors may be used
#    to endorse or promote products derived from this software without specific
#    prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
# THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
# PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL JEREMY WOHLWEND BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
# EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

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
    gene_copies = gt_params['gene_copies']
    #Remove voxel dim so that we can pass gt_params to the load_gt function
    del gt_params['voxel_dim']
    del gt_params['gene_copies']

    print "Loading data..."
    gt_dataset = load_gt(**gt_params)

    print "Labeling..."
    labeling_params = config['labeling']
    labeled_volumes, labeled_cells = label(gt_dataset, volume_dim, voxel_dim, labeling_params, gene_copies)

    print "Imaging..."
    expansion_params = config['expansion']
    optics_params = config['optics']
    volumes = resolve(labeled_volumes, volume_dim, voxel_dim, expansion_params, optics_params)
    print "Saving..."
    #Save to desired output
    output_params = config['output']
    save(volumes, **output_params)
    save_gt(gt_dataset, labeled_cells, volume_dim, volumes[0].shape, voxel_dim,\
            expansion_params, optics_params, **output_params)
    print "Done!"

if __name__ == "__main__":
    n = len(sys.argv)
    if len(sys.argv) < 2:
        raise IOError, 'Please give a config file as argument'
    #Read config file
    for i in range(1,n):
        print sys.argv[i]
        config_file = sys.argv[i]
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