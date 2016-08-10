# SimExm
A set of tools to simulate expansion microscopy experiments under different labeling and imaging conditions. 

The software is written in python with a few Cython extensions because why not ? Cython extensions require to be compiled. After installing the python dependencies you will compile the extensions use a make command as explained below. The following section covers basic usage of SimExm.


##Installation

#### Dependencies

The python dependencies can be found in the requirements.txt file. To install all at once, from the terminal, do:

`sudo pip install -r requirements.txt`  

I strongly suggest creating a virtual environment, with virtualenv or anaconda and then use the pip command to install all the dependencies inside the environment.

#### Compiling

Once all the dependencies have been correctly installed, run:  

`make install` 

which will compile the extensions. Make sure to have a look at the Makefile for useful commands.

##Basic Usage

The simulation is run through simple scripts. There are four steps to a SimExm script:  

1. Choosing a ground truth dataset to simulate.
2. Choosing labeling parameters and running the labeling simulation.
2. Choosing optics parameters and running the imaging simulation.
3. Saving, or visualizing the simulated output.

An example script can be found here: ./scripts/run_sim_easy.py
A more advanced script can be found here: ./scripts/run_sim_advanced.py

Make sure to check out both to get a better understanding of how to build a sim script.

Alternatively you can use the the IPython Notebooks which can be found in ./notebooks. 
To start the notebooks:

`cd ./notebooks; jupyter notebook`



