# SimExm
A set of tools to simulate expansion microscopy experiments under different labeling and imaging conditions. 

The software is written in python with a few Cython extensions because why not ? Cython extensions require to be compiled. After installing the python dependencies you will compile the extensions use a make command as explained below. The following section covers basic usage of SimExm.


##Installation

#### Dependencies

## Python 

I strongly suggest creating a virtual environment, with virtualenv or anaconda and then use the pip command above to install all the dependencies inside the environment. For instance, you can create a new evnironment "sim" by running:

`conda create -n sim`

Which will modify the PYTHONPATH so that all futute programs installed through pip end up in the same virtual environment.
The python dependencies can be found in the requirements.txt file. To install all at once, from the terminal, do:

`sudo pip install -r requirements.txt`  

or if from a virtual end:

\'pip install -r requirements.txt'

##Data

The data can be downloaded from : emdata.janelia.org. The following command, downloads labeled images from slice 2000 to 4000:

`for i  in {2000..4000}
do 
    wget http://emdata.janelia.org/api/node/bf1/bodies/raw/xy/2000_2000/1800_2300_$i -O destination_path/ground_truth/bodies-xy-$i.png
done`

Make sure to indicate an accessible destination path for the images

#### Compiling

Once all the dependencies have been correctly installed, run:  

`make install` or `make clean; make build`

which will cclean out the project, ompile the cython cytextensions (labeling, optics and expansion units).
Make sure to have a look at the Makefile for useful commands.

##Basic Usage

The simulation is run through simple scripts. There are four steps to a SimExm script:  

1. Creating the ground_truth object, the labeling, expansion and optics units.
2. Choosing labeling parameters and running the labeling simulation.
2. Running the imaging simulation.
3. Saving, or visualizing the simulated output.

An example script can be found here: ./scripts/run_sim.py
Make sure to check it out to get a better understanding of how to build a sim script.
The file run_sim,py can be quickly called using make run.



