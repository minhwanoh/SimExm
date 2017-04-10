channelCounts = [3, 4, 5, 6]
cellCounts = [5, 9, 13]
baseline_noises = [0, 2000, 5000]
protein_noises = [0, 0.1, 0.2]
protein_densities = [1e-6, 1e-7, 1e-8]

expansion_factor = 1.0
cellCount = cellCounts[0]
baseline_noise = baseline_noises[0]
protein_noise = protein_noises[0]
protein_density = protein_densities[0]

gene_copies = 20
fluorophores = ['ATTO488', 'ATTO550', 'ATTO700', 'Alexa790', 'ATTO490LS', 'ATTO647N']
laser_wavelength = [488, 550, 700, 790, 490, 650]
laser_filter = [[475, 550], [550, 630], [630, 720], [720, 900], [450, 550], [600, 700]]



def writeConfigFiles(cellCount, channelCount, baseline_noise, protein_noise, protein_density, expansion_factor, 
	gene_copies, fluorophores, laser_wavelength, laser_filter):
	path = "input/set3/%scells_%sch" % (cellCount, channelCount)

	groundtruth = """#Ground truth parameters
[groundtruth]

image_path = %s
offset = 49, 0, 0
bounds = 100, 200, 200
format = 'tiff'
gt_cells = 'splitted'
isotropic = False
voxel_dim = 500, 400, 400
gene_copies = %s

	[[regions]]
	#No additional regions

		""" % (path, gene_copies)

	layersHeader = """#Labeling parameters
		[labeling]
		"""

	layers = []
	for i in range(channelCount):
		layer = """
	[[layer%s]]
	fluorophore = %s
	region = 'membrane'
	labeling_density = 0.1
	protein_density = %s
	protein_noise = %s
	antibody_amp = 5.0
	single_neuron = False
		""" % (i+1, fluorophores[i], protein_density, protein_noise)
		layers.append(layer)

	expansion = """
#Expansion parameters
[expansion]

factor = %s
	""" % (expansion_factor)


	optics = """
#Optics parameters
[optics]
	
type = 'confocal'
numerical_aperture = 1.15
refractory_index = 1.33
focal_plane_depth = 500
objective_back_aperture = 1.0
exposure_time = 0.1
objective_efficiency = 0.8
detector_efficiency = 0.6
objective_factor = 40.0
pixel_size = 13500
pinhole_radius = 50.0
baseline_noise = %s

	[[channels]]
	""" % (baseline_noise)


	channels = []
	for i in range(channelCount):
		channel = """
		[[[channel_%s]]]
		laser_wavelength = %s
		laser_power = 50.0
		laser_percentage = 0.25
		laser_filter = %s
		""" % (i+1, laser_wavelength[i], str(laser_filter[i])[1:-1])
		channels.append(channel)


	output = """
[output]

name = "%scells_%sch_%sbn_%spn_%spd_%sef"
path = "output"
format =  "tiff"
sim_channels = "merged"
gt_cells = "splitted"
gt_region = "membrane"
	""" % (cellCount, channelCount, baseline_noise, protein_noise, protein_density, expansion_factor)

	filename = "config/bb_%scells_%sch_%sbn_%spn_%spd_%sef.ini" % (cellCount, channelCount, baseline_noise, protein_noise, protein_density, expansion_factor)

	# Write config file
	file = open(filename, "w") 
	file.write(groundtruth) 
	file.write(layersHeader)
	for i in range(channelCount):
		file.write(layers[i])
	file.write(expansion)
	file.write(optics)
	for i in range(channelCount):
		file.write(channels[i])
	file.write(output) 
	 
	file.close() 


expansion_factor = 1.0
baseline_noise = baseline_noises[0]
protein_noise = protein_noises[0]
protein_density = protein_densities[0]

for ch in range(len(channelCounts)):
	channelCount = channelCounts[ch]
	for cC in range(len(cellCounts)):
		cellCount = cellCounts[cC]
		writeConfigFiles(cellCount, channelCount, baseline_noise, protein_noise, protein_density, expansion_factor,
			gene_copies, fluorophores, laser_wavelength, laser_filter)


cellCount = cellCounts[1]

for ch in range(len(channelCounts)):
	channelCount = channelCounts[ch]
	for bn in range(1,len(baseline_noises)):
		baseline_noise = baseline_noises[bn]
		writeConfigFiles(cellCount, channelCount, baseline_noise, protein_noise, protein_density, expansion_factor,
			gene_copies, fluorophores, laser_wavelength, laser_filter)



baseline_noise = baseline_noises[0]

for ch in range(len(channelCounts)):
	channelCount = channelCounts[ch]
	for pn in range(1,len(protein_noises)):
		protein_noise = protein_noises[pn]
		writeConfigFiles(cellCount, channelCount, baseline_noise, protein_noise, protein_density, expansion_factor,
			gene_copies, fluorophores, laser_wavelength, laser_filter)


protein_noise = protein_noises[0]

for ch in range(len(channelCounts)):
	channelCount = channelCounts[ch]
	for pd in range(1,len(protein_densities)):
		protein_density = protein_densities[pd]
		writeConfigFiles(cellCount, channelCount, baseline_noise, protein_noise, protein_density, expansion_factor,
			gene_copies, fluorophores, laser_wavelength, laser_filter)