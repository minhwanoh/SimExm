# SimExm
A set of tools to simulate expansion microscopy experiments under different labeling and imaging conditions. 

The software is fully written in Python. The following section covers basic installation and usage of SimExm.


##Installation

The code works in python2 exclusivly. It should be easy to convert to python3 but it's untested.
I strongly suggest creating a virtual environment, with virtualenv or anaconda and then use the pip command above to install all the dependencies inside the environment. For instance, you can create a new environment "sim" by running:

`conda create -n sim python=2`

Which will modify the PYTHONPATH so that all futute programs installed through pip end up in the same virtual environment.
To start your new environment and activate the change of path run: 

`source activate sim` 

All dependencies currently used are: numpy, matplotlib, Pillow, scipy, images2gif, configobj and tifffile. 
To install all at once, from the terminal, do:

`pip install -r requirements.txt`  

or if you didn't  use a virtual environment:

`sudo pip install -r requirements.txt`

## Config Specs

The simulation is used through configuration files. We use the configobj module which provides a nice syntax and validation of the config files.
The specifications are outlines in the configspecs.ini file. Example configurations can be found under ./examples. To create a new configuration file copy one of the templates in ./examples, ad modify to fit your needs. The config parameters are explained below. For more information see the in-code documentation.

#### Groundtruth

*image_path* (string): indicates the location of the ground truth data. The data should be in the form of an image sequence, where each image represents a slice of the ground truth volume. In each image, pixel values should indicate what cell the pixel belongs to, or 0 if the pixel
is located in extra-cellular space  
*offset* (z, x, y tuple): tuple of length 3 representing where to start loading the data, in pixels. For instance, an offset of (10, 50, 50) means that images are loaded starting the 10th image in the given directory and each imaage is cropped so that the top left corner is at (50, 50).

*bounds* (z, x, y tuple): tuple of length 3 representing the size of the data to load. Could be smaller than the actual size of the ground truth. This is the size in number of voxels.
*voxel_dim*  (z, x, y tuple): tuple of length 3 representing the dimensions of a single voxel, in nanometers

*regions*: a subsection in the configuration file which may contian many different regions. A region has a single parameter, region_path, which is a string pointing to where the data for that region is. By default the software automatically computes the cytosol and membrane regions, but additional annotations may be available. They are loadede by overlapping the region with the orgiinal cell segmentaiton to figure out which part of the cell is in the given region, for each cell. See synapse.ini for an example.

#### Labeling

The labeling section should contain a subsection for each fluorophore used. See the example in brianbow_membrane.ini.
Each subsection shoudl contain the following parameters:

*fluorophore* (string) : the fluorophore to use, can be one of "ATTO390", "ATTO425", " ATTO430LS", "ATTO488", "ATTO490LS", "ATTO550", "ATTO647N", "ATTO700", "Alexa350" and "Alexa790". More info in src/fluors.py

*region* (string) : the region to annotate with the above fluorophore. Should be one of "cytosol", "membrane" or any additional regions specified in the ground truth.

*labeling_density* (float between 0.0 and 1.0): the proportion of cells in the volume to annotate with the above fluorophore

*protein_density* (float between 0.0 and 1.0): the density of proteins to label the sample with. 1.0 is 1 fluorophore per nm^3, which is not really realistic and should provide a good upperbound.

*protein_noise* (float between 0.0 and 1.0): the proportion of proteins that fly away from the labeled region. Uses a gaussian distirbution around the original location. The standard devitation for that gaussian can be modified directly in the noise function in src/labeling.py but the default should provide fairly realistic results.

*antibody_amp* (float greater than 1.0): amplifies the labeling, as well as the noise. Tipically 5.0 or 10.0.

*single_neuron* (boolean): if True, a random cell is chosen and is the only one labeled with the above fluorophore

#### Expansion

*factor* (float greater than 1.0): the expansion factor use in the expansion microscopy protocol

#### Optics
	
*numerical_aperture* (float greater than 0): the numerical aperture of the system
*focal_plane_depth* (integer greater than 1): the thickness of a z-slice in nanometers
*objective_back_aperture* (float greater than 0): the objective back aperture of the system
*exposure_time* (float greater than 0): how long photons are detected in seconds
*objective_efficiency* (float between 0.0 and 1.0): the percentage of efficiency of the objective
*detector_efficiency* (float between 0.0 and 1.0): the percentage of efficiency of the detector
*objective_factor* (float greater or equal to 0.0): the objective lens, tipically 20x, 40x
*pixel_size* (integer): the size of an output pixel in the microscope, in nanometers
*baseline_noise* (integer): the average number of baseline photons detected by the system
	
*channels* : subsection containing a multiple channel parameters for different lasers. Each subsection has the following parameters. See brainbow_membrane.ini for an example on how to use multiple channels.

*laser_wavelength* (ineteger between 200 and 1000): the wavelength of the laser in nanometers
*laser_power* (float greater than 1.0): the power of the laser in Watts
*laser_percentage* (float between 0.0 and 1.0): proportion of the laser power to use
*laser_filter* (min, max integer tuple): the minimum and maximum wavelenghts in the filter, in nanometers

#### Output

The simualtion outputs a simualted stack and the corresponding ground truth, which can be used for validation or training.
Output formats include: tiff stacks, gif or png sequence. The simulation stack may be outputted in separate channels or merged
into multiple RGB volumes. If the number of channels is greater than 3, an RGB volume may be created for every 3 channels.

The ground truth is stored on a per fluorophore basis. For each fluorophore the output may be spliited into a volume for each cell
or grouped into a single volume containing all cells labeled by that fluorophore. The parameters are outlined below.

*name* (string): a name for the experiment
*path* (path): where to store the experiement's output
*format* (string): may be one of 'tiff', 'gif' or 'image sequence', determines to output format for both the simulated stack and the ground truth
*sim_channels* (string): may be one of "merged" or "splitted", if merged, and RGB volume is made for every 3 channels, otherwise a stack is made for each channel
*gt_cells* (string): may be one of "merged ot "splitted", if merged, all cells are grouped in the same stack, otherwise, a volume if made for each cell
*gt_region* (string): the cell region to use in the output, may be different that the annotated regions in the labeling layers.

## Run the Simulation

Once the config.ini file is ready, make sure that you have activated your virtual environment if you have one (using source activate environment_name). Then, run the simulation from the terminal by using the following command:

`python run.py path_to_config/config.ini`

Note that the simulation may require a decent amount of memory for large volumes. A computer with 4G of RAM can usually handle a volume of about 500 x 500 x 500 voxels without problems. If you run into memory issues, consider using a computer with more RAM.





