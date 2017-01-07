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
optics.py

Set of methods to handle the simulation of the optical process.
Implements confocal light microscopy.
"""

import numpy as np
from psf import gaussian_psf, wolf_born_psf
from fluors import Fluorset
from scipy.signal import fftconvolve
from scipy.misc import imresize

def resolve(labeled_volumes, volume_dim, voxel_dim, expansion_params, optics_params):
    """
    Resolves the labeled volumes with the given optics parameters.
    Performs photon count calculation, convolution with a point spread function,
    baseline noise and rescaling.

    Args:
        labeled_volumes: dict fluorophore (string) -> volume (numpy 3D uint32 array)
            dictionary containing the volumes to resolve
        volume_dim: (z, x, y) integer tuple
            dimensions of each volume in number of voxels
        voxel_dim: (z, x, y) integer tuple
            dimensions of a voxel in nm
        expansion_parameters: dict
            dicitonary containing the expansion parameters
        optics_parameters: dict
            dicitonary containing the optics parameters
    Returns:
        volumes: list of numpy 3D uint8 arrays
            as list contianing a volume resolved for each channel
    """
    #Create volume
    volumes = []
    #Resolve each channel one by one
    #Make sure they're sorted by name for consistency
    channels = sorted(optics_params['channels'].keys())
    for channel in channels:
        print "Resolving {}".format(channel)
        vol = np.zeros(volume_dim, np.uint32)
        channel_params = optics_params['channels'][channel]
        #Compute photon count
        for fluorophore in labeled_volumes:
            params = optics_params.copy()
            params.update(channel_params)
            mean_photon = mean_photons(fluorophore, **params)
            Z, X, Y = np.nonzero(labeled_volumes[fluorophore])
            photons = np.random.poisson(mean_photon, size = len(Z)).astype(np.uint32)
            photons = np.multiply(labeled_volumes[fluorophore][Z, X, Y], photons)
            np.add.at(vol, (Z, X, Y), photons)
        #Convolve with point spread
        psf_voxel_dim = np.array(voxel_dim) * expansion_params['factor']
        psf = gaussian_psf(psf_voxel_dim, **params)
        out = np.round(fftconvolve(vol.astype(np.float64), psf, 'valid'))
        #Add noise
        out = np.add(out, baseline_volume(out.shape, **optics_params))
        #Optical scaling
        out = scale(out, voxel_dim, 'bilinear', expansion_params['factor'], **optics_params)
        #Normalize
        out = normalize(out)
        volumes.append(out)

    return volumes

def baseline_volume(volume_dim, baseline_noise, **kwargs):
    """
    Creates a volume of baseline photon noise, using a poisson distribution
    multiplied by a gaussian to mimic the light focus towards the center of the image.

    Args:
        volume_dim: (z, x, y) integer tuple
            the size of the volume to create
        baseline_noise: integer 
            the mean number of photons of the poisson distribution
    Returns:
        out: numpy 3D float64 arrat
            a volume of basline photon noise
    """
    (d, w, h) = volume_dim
    indices = np.indices((d, w, h), dtype=np.float64)
    gaussian = np.exp(-((indices[1] - w / 2.0)**2 / (0.5 * w**2) + (indices[2] - h / 2.0)**2 / (0.5 * h**2)))
    out = np.round(np.multiply(np.random.poisson(baseline_noise, size=volume_dim), gaussian))
    return out

def normalize(volume):
    """
    Normalizes the volume to the [0, 255] range by dividing by the maximum value.
    
    Args:
        volume: 3D numpy array, np.uint32
            the volume to normalize
    Returns:
        normalized:3D numpy array, np.uint8
            the normalized volume
    """
    max = np.amax(volume)
    if max == 0:#Fixes dividing by 0 error if nothing in the volume
        return volume.astype(np.uint8)
    normalized = volume * (255.0 / max)
    normalized = np.round(normalized).astype(np.uint8)
    return normalized

def scale(volume, voxel_dim, interpolation, expansion, objective_factor,
            pixel_size, focal_plane_depth, **kwargs):
    """
    Scales the output volume with the appropriate optics parameters
    using nearest neighbour interpolation.
    
    Args:
        volume: numpy 3D array (z, x, y)
            the volume to scale
        voxel_dim: (z, x, y) tuple
            the dimensions of a voxel in nm
        interpolation: string
            one of 'nearest', 'bilinear', 'bicubic' and 'cubic'
        expansion: float
            the expansion factor
        objective_factor: float
            objective factor of the microscope, tipically 0, 20 or 40
        pixel_size: integer
            the size of a pixel in the microscope, in nm
        focal_plane_depth: integer
            the thickness of a slice in nm
    Returns:
        out: numpt 3D array
            the scaled volume in all three axis
    """
    xy_scale  = (voxel_dim[1] * expansion * objective_factor)  / float(pixel_size)
    z_scale = voxel_dim[0] * expansion / float(focal_plane_depth)
    z_step = np.ceil(1.0 / z_scale).astype(np.int)
    out = []
    for i in range(0, volume.shape[0], z_step):
        im = imresize(volume[i, : ,:], xy_scale, interp = interpolation)
        out.append(im)
    return np.array(out)

def mean_photons(fluorophore, exposure_time, objective_efficiency,\
                detector_efficiency, objective_back_aperture, objective_factor, \
                laser_wavelength, laser_filter, laser_power, laser_percentage, **kwargs):
    """
    Computes the mean number of detected photons for a given fluorophore and laser parameters

    Args:
        fluorophore: string
            the fluorophore to measure
        exposure_time: float
            amoung of time that the laser is shined on the specimen, in seconds
        objective_efficiency: float
            the efficiency of the microscope's objective
        detector_efficiency: float
            the efficiency of the microsocpe's detector
        objective_back_aperture: float
            name says it all
        objective_factor:
            objective lense, tipically 0, 20 or 40
        laser_wavelength: integer
            the wavelength of the laser, in nm
        laser_filter: integer list of length 2 
            wavelength_min and wavelength_maximum detected
        laser_power: float
            the laser power
        laser_percentage: float
            proportion of power to use
    Returns:
        detected_photons: integer
            the mean number of detected photons per fluorophore protein
    """
    fluorset = Fluorset()
    #Get fluorophore data
    f = fluorset.get_fluor(fluorophore)
    qy, ext_coeff = f.get_quantum_yield(), f.get_extinction_coefficient()
    excitation, emission = f.find_excitation(laser_wavelength), f.find_emission(laser_filter)
    #Compute parameters
    laser_radius = float(objective_back_aperture) / objective_factor
    laser_intensity = laser_power * laser_percentage / (np.pi * laser_radius**2)# [Watts / (m ^ 2)]
    CONSTANT = 0.119626566 # in m^3 * kg * s^{-1} Avogadro's number * Planck's constant * speed of light
    #Get mean photon count
    emitted_photons = excitation * qy * (ext_coeff * 1e2) * exposure_time *\
                      laser_intensity * (laser_wavelength * 1e-9) / (1e3 * CONSTANT)
    detected_photons = emitted_photons * emission * objective_efficiency * detector_efficiency
    detected_photons = int(np.round(detected_photons))
    return detected_photons
