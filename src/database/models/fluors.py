'''
fluors.py

The Fluorophore class and a subclass for each type (for easy data access and storage)

'''

import os
import numpy as np


class Fluorophore:

    '''
    Interface to interace with a fluorophore's data.
    Instances of the FLuorohore class are used to enumerate fluorophore parameters,
    some of which are hard-coded in this file.

    Init method, sets attributes

    name (string) : Name given to the fluorophore
    xls (string) : Path to xls file containing fluorophore data
    excitation (string) : string, Path to file containing excitation data
    emission (string) : Path to file containing excitation data
    extinction_coefficient (int) : Coefficient of exctinction
    quantum_yield (float) : Value of the quantum yeild
    source (string) : Source where the data can be obtained / verified
    comments (string) : Additonal comments for the user 
    '''

    name = ""
    xls = ""
    quantum_yield = 0.0
    extinction_coefficient = 0.0
    emission = ""
    excitation = ""
    source = ""
    comments = ""

    def get_name(self):
        ''' Returns the (string) type of the fluoropore '''

        return self.name

    def get_quantum_yield(self):
        ''' Returns the (float) quantum yield of the fluorophore '''

        return float(self.quantum_yield)

    def get_extinction_coefficient(self):
        ''' Returns the (float) extinction coefficient of the fluorophore '''

        return float(self.extinction_coefficient)

    def get_source(self):
        '''Returns the (string) source, where information about the flurophore can be found '''

        return self.source

    def get_comments(self):
        ''' Returns (string) comments, addional information about the fluorophore '''

        return self.comments

    def get_emission_file(self):
        ''' Returns emission data as a float numpy array '''

        file_path = self.emission
        f = open(file_path, 'r')
        raw_data = f.read().split("\r")
        raw_data = np.array([s.split("\t") for s in raw_data], dtype=float)
        return raw_data

    def get_excitation_file(self):
        ''' Returns excitation data as a float numpy array '''

        file_path = self.excitation
        f = open(file_path, 'r')
        raw_data = f.read().split("\r")
        raw_data = np.array([s.split("\t") for s in raw_data], dtype=float)
        return raw_data

    def find_excitation(self, wavelength):
        ''' 
        Returns the excitation value of the given fluorophore, at the given wavelenght

        wavelenght (int) : value in nm
        '''

        raw_data = self.get_excitation_file()
        data     = raw_data[:,1] / np.max(raw_data[:,1])
        data_min = np.min(raw_data[:,0])
        data_max = np.max(raw_data[:,0])
        interval = round(np.mean(np.diff(raw_data[:,0])),2)

        if (wavelength < data_min) or (wavelength > data_max):
            excitation = 0
        else:
            index      = int(np.floor((wavelength - data_min) / interval))
            index -= 1
            weight     = wavelength - (data_min + index * interval)
            excitation = weight * data[index + 2] + (1 - weight) * data[index + 1]

        return excitation

    def find_emission(self, wavelength_min, wavelength_max):
        ''' 
        Returns the emission value of the given fluorophore, at the given wavelenghts

        wavelenght_min (int) : minimum value in nm
        wavelenght_max (int) : maximum value in nm
        '''

        raw_data = self.get_emission_file()
        data     = raw_data[:,1] / np.sum(raw_data[:,1])
        interval = round(np.mean(np.diff(raw_data[:,0])),2)
        data_min = np.min(raw_data[:,0])
        data_max = np.max(raw_data[:,0])
        emission = 0

        if wavelength_max < data_max:
            index = int(np.floor((wavelength_max - data_min) / interval))
        else:
            index = len(data) - 1

        index -= 1

        while (index >= 0) and (wavelength_min <= data_min + index * interval):
            emission = emission + data[index + 1]
            index = index - 1

        return emission

#Path to data
fluo_path = os.path.abspath(os.path.dirname(__file__)) + '/../fluordata/'

# List of fluorophores in database, as FLuorophore objects
# Alexa DYES

class Alexa350(Fluorophore):

    name = "Alexa350"
    xls = fluo_path  + "Alexa350.xls"
    excitation = fluo_path +"Alexa350_excitation.txt"
    emission = fluo_path +"Alexa350_emission.txt"
    extinction_coefficient = 19000
    quantum_yield = 0.25 
    source = "https://www.thermofisher.com/us/en/home/references/molecular-probes-the-handbook/fluorophores-and-their-amine-reactive-derivatives/alexa-fluor-dyes-spanning-the-visible-and-infrared-spectrum.html"
    comments = "The quantum yeild value is a guess, the actual figure is not provided."

class Alexa790(Fluorophore):

    name = "Alexa790"
    xls = fluo_path  + "Alexa790.xls"
    excitation = fluo_path +"Alexa790_excitation.txt"
    emission = fluo_path +"Alexa790_emission.txt"
    extinction_coefficient = 260000
    quantum_yield = 0.25 
    source = "https://www.thermofisher.com/us/en/home/references/molecular-probes-the-handbook/fluorophores-and-their-amine-reactive-derivatives/alexa-fluor-dyes-spanning-the-visible-and-infrared-spectrum.html"
    comments = "The quantum yeild value is a guess, the actual figure is not provided."

# ATTO DYES

class ATTO390(Fluorophore):

    name = "ATTO390"
    xls = fluo_path  + "ATTO390.xls"
    excitation = fluo_path +"ATTO390_excitation.txt"
    emission = fluo_path +"ATTO390_emission.txt"
    extinction_coefficient = 24000
    quantum_yield = 0.9 
    source = "http://www.atto-tec.com/fileadmin/user_upload/Katalog_Flyer_Support/ATTO_390.pdf"
    comments = ""

class ATTO425(Fluorophore):

    name = "ATTO425"
    xls = fluo_path  + "ATTO425.xls"
    excitation = fluo_path +"ATTO425_excitation.txt"
    emission = fluo_path +"ATTO425_emission.txt"
    extinction_coefficient = 45000
    quantum_yield = 0.9 
    source = "http://www.atto-tec.com/fileadmin/user_upload/Produktdatenblaetter/ATTO_425.pdf"
    comments = ""

class ATTO430LS(Fluorophore):

    name = "ATTO430LS"
    xls = fluo_path  + "ATTO430LS.xls"
    excitation = fluo_path +"ATTO430LS_excitation.txt"
    emission = fluo_path +"ATTO430LS_emission.txt"
    extinction_coefficient = 32000
    quantum_yield = 0.65 
    source = "http://www.atto-tec.com/fileadmin/user_upload/Katalog_Flyer_Support/ATTO_430LS.pdf"
    comments = ""

class ATTO465(Fluorophore):

    name = "ATTO465"
    xls = fluo_path  + "ATTO465.xls"
    excitation = fluo_path +"ATTO465_excitation.txt"
    emission = fluo_path +"ATTO465_emission.txt"
    extinction_coefficient = 75000
    quantum_yield = 0.75 
    source = "http://www.atto-tec.com/fileadmin/user_upload/Katalog_Flyer_Support/ATTO_465.pdf"
    comments = ""

class ATTO488(Fluorophore):

    name = "ATTO488"
    xls = fluo_path  + "ATTO488.xls"
    excitation = fluo_path +"ATTO488_excitation.txt"
    emission = fluo_path +"ATTO488_emission.txt"
    extinction_coefficient = 90000
    quantum_yield = 0.8 
    source = "http://www.atto-tec.com/fileadmin/user_upload/Katalog_Flyer_Support/ATTO_488.pdf"
    comments = ""

class ATTO490LS(Fluorophore):

    name = "ATTO490LS"
    xls = fluo_path  + "ATTO490LS.xls"
    excitation = fluo_path +"ATTO490LS_excitation.txt"
    emission = fluo_path +"ATTO490LS_emission.txt"
    extinction_coefficient = 40000
    quantum_yield = 0.3 
    source = "http://www.atto-tec.com/fileadmin/user_upload/Katalog_Flyer_Support/ATTO_490LS.pdf"
    comments = ""

class ATTO550(Fluorophore):

    name = "ATTO550"
    xls = fluo_path  + "ATTO550.xls"
    excitation = fluo_path +"ATTO550_excitation.txt"
    emission = fluo_path +"ATTO550_emission.txt"
    extinction_coefficient = 120000
    quantum_yield = 0.3 
    source = "http://www.atto-tec.com/fileadmin/user_upload/Katalog_Flyer_Support/ATTO_550.pdf"
    comments = ""

class ATTO647N(Fluorophore):

    name = "ATTO647N"
    xls = fluo_path  + "ATTO647N.xls"
    excitation = fluo_path +"ATTO647N_excitation.txt"
    emission = fluo_path +"ATTO647N_emission.txt"
    extinction_coefficient = 150000
    quantum_yield = 0.65 
    source = "http://www.atto-tec.com/fileadmin/user_upload/Katalog_Flyer_Support/ATTO_647N.pdf"
    comments = ""

class ATTO700(Fluorophore):

    name = "ATTO700"
    xls = fluo_path  + "ATTO700.xls"
    excitation = fluo_path +"ATTO700_excitation.txt"
    emission = fluo_path +"ATTO700_emission.txt"
    extinction_coefficient = 120000
    quantum_yield = 0.25 
    source = "http://www.atto-tec.com/fileadmin/user_upload/Katalog_Flyer_Support/ATTO_700.pdf"
    comments = ""
