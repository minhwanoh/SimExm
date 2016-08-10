import os

"""
List of fluorophores avalable for labeling simulation.
Each fluorophore has the following parameters:

@name (string) : Name given to the fluorophore
@xls (string) : Path to xls file containing fluorophore data
@excitation (string) : string, Path to file containing excitation data
@emission (string) : Path to file containing excitation data
@extinction_coefficient (int) : Coefficient of exctinction
@quantum_yield (float) : Value of the quantum yeild
@source (string) : Source where the data can be obtained / verified
@comments (string) : Additonal comments for the user 

"""

fluo_path = os.path.abspath(os.path.dirname(__file__)) + '/FluorData/'

# Alexa DYES

class Alexa350:

    name = "Alexa350"
    xls = fluo_path  + "Alexa350.xls"
    excitation = fluo_path +"Alexa350_excitation.txt"
    emission = fluo_path +"Alexa350_emission.txt"
    extinction_coefficient = 19000
    quantum_yield = 0.25 
    source = "https://www.thermofisher.com/us/en/home/references/molecular-probes-the-handbook/fluorophores-and-their-amine-reactive-derivatives/alexa-fluor-dyes-spanning-the-visible-and-infrared-spectrum.html"
    comments = "The quantum yeild value is a guess, the actual figure is not provided."

class Alexa790:

    name = "Alexa790"
    xls = fluo_path  + "Alexa790.xls"
    excitation = fluo_path +"Alexa790_excitation.txt"
    emission = fluo_path +"Alexa790_emission.txt"
    extinction_coefficient = 260000
    quantum_yield = 0.25 
    source = "https://www.thermofisher.com/us/en/home/references/molecular-probes-the-handbook/fluorophores-and-their-amine-reactive-derivatives/alexa-fluor-dyes-spanning-the-visible-and-infrared-spectrum.html"
    comments = "The quantum yeild value is a guess, the actual figure is not provided."

# ATTO DYES

class ATTO390:

    name = "ATTO390"
    xls = fluo_path  + "ATTO390.xls"
    excitation = fluo_path +"ATTO390_excitation.txt"
    emission = fluo_path +"ATTO390_emission.txt"
    extinction_coefficient = 24000
    quantum_yield = 0.9 
    source = "http://www.atto-tec.com/fileadmin/user_upload/Katalog_Flyer_Support/ATTO_390.pdf"
    comments = ""

class ATTO425:

    name = "ATTO425"
    xls = fluo_path  + "ATTO425.xls"
    excitation = fluo_path +"ATTO425_excitation.txt"
    emission = fluo_path +"ATTO425_emission.txt"
    extinction_coefficient = 45000
    quantum_yield = 0.9 
    source = "http://www.atto-tec.com/fileadmin/user_upload/Produktdatenblaetter/ATTO_425.pdf"
    comments = ""

class ATTO430LS:

    name = "ATTO430LS"
    xls = fluo_path  + "ATTO430LS.xls"
    excitation = fluo_path +"ATTO430LS_excitation.txt"
    emission = fluo_path +"ATTO430LS_emission.txt"
    extinction_coefficient = 32000
    quantum_yield = 0.65 
    source = "http://www.atto-tec.com/fileadmin/user_upload/Katalog_Flyer_Support/ATTO_430LS.pdf"
    comments = ""

class ATTO465:

    name = "ATTO465"
    xls = fluo_path  + "ATTO465.xls"
    excitation = fluo_path +"ATTO465_excitation.txt"
    emission = fluo_path +"ATTO465_emission.txt"
    extinction_coefficient = 75000
    quantum_yield = 0.75 
    source = "http://www.atto-tec.com/fileadmin/user_upload/Katalog_Flyer_Support/ATTO_465.pdf"
    comments = ""

class ATTO488:

    name = "ATTO488"
    xls = fluo_path  + "ATTO488.xls"
    excitation = fluo_path +"ATTO488_excitation.txt"
    emission = fluo_path +"ATTO488_emission.txt"
    extinction_coefficient = 90000
    quantum_yield = 0.8 
    source = "http://www.atto-tec.com/fileadmin/user_upload/Katalog_Flyer_Support/ATTO_488.pdf"
    comments = ""

class ATTO490LS:

    name = "ATTO490LS"
    xls = fluo_path  + "ATTO490LS.xls"
    excitation = fluo_path +"ATTO490LS_excitation.txt"
    emission = fluo_path +"ATTO490LS_emission.txt"
    extinction_coefficient = 40000
    quantum_yield = 0.3 
    source = "http://www.atto-tec.com/fileadmin/user_upload/Katalog_Flyer_Support/ATTO_490LS.pdf"
    comments = ""

class ATTO550:

    name = "ATTO550"
    xls = fluo_path  + "ATTO550.xls"
    excitation = fluo_path +"ATTO550_excitation.txt"
    emission = fluo_path +"ATTO550_emission.txt"
    extinction_coefficient = 120000
    quantum_yield = 0.3 
    source = "http://www.atto-tec.com/fileadmin/user_upload/Katalog_Flyer_Support/ATTO_550.pdf"
    comments = ""

class ATTO647N:

    name = "ATTO647N"
    xls = fluo_path  + "ATTO647N.xls"
    excitation = fluo_path +"ATTO647N_excitation.txt"
    emission = fluo_path +"ATTO647N_emission.txt"
    extinction_coefficient = 150000
    quantum_yield = 0.65 
    source = "http://www.atto-tec.com/fileadmin/user_upload/Katalog_Flyer_Support/ATTO_647N.pdf"
    comments = ""

class ATTO700:

    name = "ATTO700"
    xls = fluo_path  + "ATTO700.xls"
    excitation = fluo_path +"ATTO700_excitation.txt"
    emission = fluo_path +"ATTO700_emission.txt"
    extinction_coefficient = 120000
    quantum_yield = 0.25 
    source = "http://www.atto-tec.com/fileadmin/user_upload/Katalog_Flyer_Support/ATTO_700.pdf"
    comments = ""

#Create accessors

types = [Alexa350, Alexa790, ATTO390, ATTO425, ATTO430LS, ATTO465, ATTO488, ATTO490LS, ATTO550, ATTO647N, ATTO700]
name_map = {"Alexa350": Alexa350, "Alexa790": Alexa790, "ATTO390": ATTO390, "ATTO425": ATTO425, "ATTO430LS": ATTO430LS,\
 "ATTO465": ATTO465, "ATTO488": ATTO488, "ATTO490LS": ATTO490LS, "ATTO550": ATTO550, "ATTO647N": ATTO647N, "ATTO700": ATTO700}
