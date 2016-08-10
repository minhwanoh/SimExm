import numpy as np
from math import floor
import database.externals.fluors as fdata

############# FLUOROPHORE #################

""" Get the location of the fluorphore given its id """
def get_all_fluorophore_types(ds):
    return fdata.types

def get_quantum_yield(ds, name):
    fluor = fdata.name_map[name]
    return float(fluor.quantum_yield)

def get_extinction_coefficient(ds, name):
    fluor = fdata.name_map[name]
    return float(fluor.extinction_coefficient)

def get_emission_file(ds, fluo_name):
    fluor = fdata.name_map[fluo_name]
    file_path = fluor.emission
    f = open(file_path, 'r')
    raw_data = f.read().split("\r")
    raw_data = np.array([s.split("\t") for s in raw_data], dtype=float)
    return raw_data

def get_excitation_file(ds, fluo_name):
    fluor = fdata.name_map[fluo_name]
    file_path = fluor.excitation
    f = open(file_path, 'r')
    raw_data = f.read().split("\r")
    raw_data = np.array([s.split("\t") for s in raw_data], dtype=float)
    return raw_data

''' Returns the excitation value of the given fluorophore, at the given wavelenght '''
def find_excitation(ds, fluo_name, wavelength):
    raw_data = get_excitation_file(ds, fluo_name)
    data     = raw_data[:,1] / np.max(raw_data[:,1])
    data_min = np.min(raw_data[:,0])
    data_max = np.max(raw_data[:,0])
    interval = round(np.mean(np.diff(raw_data[:,0])),2)

    if (wavelength < data_min) or (wavelength > data_max):
        excitation = 0
    else:
        index      = int(floor((wavelength - data_min) / interval))
        index -= 1
        weight     = wavelength - (data_min + index * interval)
        excitation = weight * data[index + 2] + (1 - weight) * data[index + 1]

    return excitation

''' Returns the emission value of the given fluorophore, at the given wavelenght '''
def find_emission(ds, fluo_name, wavelength_min, wavelength_max):
    raw_data = get_emission_file(ds, fluo_name)
    data     = raw_data[:,1] / np.sum(raw_data[:,1])
    interval = round(np.mean(np.diff(raw_data[:,0])),2)
    data_min = np.min(raw_data[:,0])
    data_max = np.max(raw_data[:,0])
    emission = 0

    if wavelength_max < data_max:
        index = int(floor((wavelength_max - data_min) / interval))
    else:
        index = len(data) - 1

    index -= 1

    while (index >= 0) and (wavelength_min <= data_min + index * interval):
        emission = emission + data[index + 1]
        index = index - 1

    return emission

