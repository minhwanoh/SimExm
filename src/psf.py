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
psf.py

Contains different point spread function models for simulation purposes
"""

import numpy as np
from itertools import product
from scipy.integrate import quad as integrate
from scipy.special import j0

def gaussian_psf(voxel_dim, numerical_aperture, laser_wavelength, **kwargs):
    #e^2 - 1 = 6.389
    shape = np.ceil(800.0 / np.array(voxel_dim)).astype(np.int)
    #Compute size coefficient
    d, w, h = voxel_dim
    #Set bounds
    indices = np.indices(shape, dtype=np.float64)
    #Compute using gaussian beam formula
    w_0 = laser_wavelength / (np.pi * numerical_aperture)
    z_r = np.pi * ((w_0)**2) / laser_wavelength
    w_z = (w_0) * np.sqrt(1 + (d * indices[0] / z_r)**2)
    kernel = (w_0 / w_z) * np.exp(- ((w * indices[1])**2 + (h * indices[2])**2) / (w_z ** 2))
        
    return mirror(kernel)

def wolf_born_psf(voxel_dim, refr_index, numerical_aperture, laser_wavelength, **kwargs):
    d, w, h = voxel_dim
    w_0 = laser_wavelength / (np.pi * numerical_aperture)
    z_r = np.pi * ((w_0)**2) / laser_wavelength
    #e^2 - 1 = 6.389
    s = int(np.ceil(6.389 * z_r / d))
    psf = np.zeros((s, s, s), np.float64)
    for k, i, j in product(xrange(s), xrange(s), xrange(s)):
        def f1(rho):
            tmp = j0(numerical_aperture * 2 * np.pi/ (refr_index * laser_wavelength) * ((w*i)**2 + (h*j)**2)**0.5 * rho)
            return tmp * np.cos(- np.pi / laser_wavelength * rho**2 * d * k * (numerical_aperture / refr_index)**2)
        def f2(rho):
            tmp = j0(numerical_aperture  * 2 * np.pi / (refr_index * laser_wavelength) * ((w*i)**2 + (h*j)**2)**0.5 * rho)
            return tmp * np.sin(- np.pi / laser_wavelength * rho**2 * d * k * (numerical_aperture / refr_index)**2)
        psf[k, i, j] = integrate(f1, 0.0, 1.0)[0]**2 + integrate(f2, 0.0, 1.0)[0]**2 
    return mirror(psf)

def mirror(vol):
    out = np.zeros([2*i-1 for i in vol.shape], np.float64)
    x, y, z = (i - 1 for i in vol.shape)
    out[x:, y:, z:] = vol
    out[:x, y:, z:] = vol[-1:0:-1, :, :]
    out[:, :y, z:] = out[:, -1:y:-1, z:]
    out[:, :, :z] = out[:, :, -1:z:-1]
    return out

def normalize(vol):
    return 1.0 / np.sum(vol) * vol