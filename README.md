# SimExm
A set of tools to simulate expansion microscopy experiments under different labeling and imaging conditions. 

The software is written in python with a few Cython extensions because why not ? Cython extensions require to be compiled. After installing the python dependencies you will compile the extensions using a make command as explained below. The following section covers basic installation and usage of SimExm.


##Installation

#### Dependencies

The code works in python2 exclusivly. It should be easy to convert to python3 but it's untested.
I strongly suggest creating a virtual environment, with virtualenv or anaconda and then use the pip command above to install all the dependencies inside the environment. For instance, you can create a new environment "sim" by running:

`conda create -n sim python`

Which will modify the PYTHONPATH so that all futute programs installed through pip end up in the same virtual environment.
To start your new environment and activate the change of path run: 

`source activate sim` 

All dependencies currently used are: numpy, Cython, jupyter, ipywidgets, matplotlib, Pillow, scipy, images2gif,and h5py. 
To install all at once, from the terminal, do:

`pip install -r requirements.txt`  

or if you didn't  use a virtual environment:

`sudo pip install -r requirements.txt`

#### Dataset

I'm working on getting better ground truth data. In the mean time I use a dataset from the Janelia group, which unfortuntaly seems to contain erros of annotation. It is, however, good enough for the purpose of these simulations.
The data can be downloaded from : emdata.janelia.org. The following command, downloads labeled images from slice 2000 to 4000:

```
for i  in {2000..4000}  
do  
    wget http://emdata.janelia.org/api/node/bf1/bodies/raw/xy/2000_2000/1800_2300_$i -o destination_path/ground_truth/bodies-xy-$i.png  
done
```

Make sure to indicate a writable destination path for the images.

#### Compiling

Once all the dependencies have been correctly installed, run:  

`make install` or `make clean; make build`

which will clean out the project, compile the cython extensions (labeling, optics and expansion units).
Make sure to have a look at the Makefile for useful commands.

##Basic Usage

The simulation is run through simple scripts. There are 6 steps to a SimExm script:  

1. Creating the ground_truth dataset object, as well as the labeling, 
2. Create the expansion and optics units and specify parameters.
3. Fill up the ground truth dataset object with data from a saved file or by loading from the images.
4. Choosing labeling parameters and running the labeling simulation.
5. Running the expansion and imaging simulation.
6. Saving, or visualizing the simulated output.

An example script can be found here: ./examples/run_sim.py
Make sure to check it out to get a better understanding of how to build a sim script.
The file run_sim,py can be quickly called using `make run`.



