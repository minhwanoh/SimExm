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
output.py

Handles storage of simulation outputs.
Simulation stacks can be saved in three possible formats:
	- TIFF stack
	- GIF stack
	- Image sequence (png)
In addition, one can also merge channels into rgb volumes or save each channel
separatly. The same can be done for the ground truth output
"""

import os
from tifffile import imsave

SAVE_FUNCTION = {'tiff': save_as_tiff,\
		  		 'gif': save_as_gif,\
	      		 'image sequence': save_as_image_sequence}

def save(volumes, path, name, sim_channels, format, **kwargs):
	"""
	Saves the simulation stack with the given output parameters.

	Args:
		volumes
		path
		name
		sim_channels
		format
	Returns:
		None
	"""
	if path[-1] != "/": path.append("/")
	os.mkdir(path + name)

	if sim_channels == 'merged':
		volumes = merge(volumes)
		i = 0
		for vol in volumes:
			sf = SAVE_FUNCTION[format]
			sf(vol, path=path + name + '/', 'channels_{}{}{}'.format(i, i + 1, i + 2))
			i += 3
	else:
		i = 0
		for vol in volumes:
			sf = SAVE_FUNCTION[format]
			sf(vol, path=path + name + '/', 'channel_{}'.format(i))
			i += 1

def save_gt(labeled_cells, scaling_factor, path, name, gt_channels, format, **kwargs):
	"""
	Saves the ground truth stack with the given output parameters
	"""

def merge(volumes):
	"""
	Merges the given volumes into multiple 3 channel volumes

	Args:
		volumes: list of numpy 3D arrays
			list of volumes to break up and stack into multiple 3-channels stacks
	Returns:
		out: list of numpy 4D arrays (z, x, y, channel)
			list of rgb volumes
	"""
	out = []
	#Add missing channels to round up to % 3
	for i in xrange(len(volumes) % 3)
		(d, w, h) = volumes[0].shape
		empty = np.zeros_like(volumes[0])
		out.append(empty)
	#Split every 3 stacks
	for i in xrange(0, len(volumes), 3):
		vol = np.stack(volumes[i:i+3], axis = -1)
		out.append(vol)
	return out

def save_as_tiff(volume, path, name, rgb):
	dest = path + name + '.tiff'
	if rbg:
		imsave(dest, volume, photometric='rgb')
	else:
		imsave(dest, volume, photometric='minisblack')

def save_as_gif(volume, path, name, rgb):

def save_as_image_sequence(volume, path, name, rgb):

