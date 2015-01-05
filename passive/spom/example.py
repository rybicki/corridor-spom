"""
 An example script for running SPOM simulations.
"""

import time
import numpy

import model.grid
from model.grid import Species, Grid
import habitat

# The size of the simulation grid
height, width = 256, 256

# Simulation time
time_steps = 100

# A flag to indicate whether to use matplotlib the results. 
# If you do not have matplotlib, set this to False.
do_plot = True

# Initialize the C module 
print "Initializing the SPOM library"
threads = 1 # number of threads used by FFTW library
seed = int(time.time()) # the random seed to be used in the simulation
model.grid.initialize_model(threads, seed)

# Create a random autocorrelated habitat matrix. 
print "Generating a habitat matrix"
habitat_matrix = habitat.generate_habitat(height, width, 1.5)

# A list of species to be simulated. Species parameters are given in the following order:
# - phenotype (phi)
# - colonization rate (c)
# - extinction rate (e)
# - inverse of dispersal distance (alpha)
# - niche width (gamma)
species_list = [
  Species(0.8, 0.5, 0.1, 0.1, 0.1),
  Species(0.1, 2.0, 0.3, 1.0, 0.2),
  Species(0.5, 0.75, 0.2, 0.5, 0.1)
]

def make_spom_grid(species):
	# Each species needs it own grid instance
	g = Grid(height, width, species)
	g.initialize()

	# All species share the same habitat matrix
	g.calculate_fitness_from_habitat(habitat_matrix) 
	
	# Initialize the occupancy matrix to contain all 1s.
	# That is, all patches are initially occupied.
	g.occupancy[:] = 1 	

	return g
	
# Initialize a new SPOM for each species. 
print "Initializing SPOMs for the specified species parameters"
spoms = [make_spom_grid(sp) for sp in species_list]

# Run the simulation for each grid
histories = [] 
for (index, spom) in enumerate(spoms):
	print "Running simulation %i for %i steps or until extinction.." % (index, time_steps,)
	begin_time = time.clock()
	history = spom.many_steps(time_steps) # Returns number of occupied patches during each step
	end_time = time.clock()
	print "Simulation took %.3f seconds" % (end_time - begin_time)
	histories.append(history)

# ---------- Plot the simulation results ----------
# This section requires the 'pylab' module that comes with 'matplotlib' library.
# In case you want to output the simulation results some other way, you do not need to
# install the 'matplotlib' library.

def plot_habitat(habitat_fig):
	# Show the random habitat matrix
	pylab.title("Habitat matrix")
	pylab.imshow(habitat_matrix)
	print "Writing '%s'" % (habitat_fig,)
	pylab.savefig(habitat_fig)
	pylab.show()

def plot_occupancy(spom, occupancy_fig):
	# Plot the final occupancy matrix
	pylab.title("Occupancy matrix")
	pylab.imshow(spom.occupancy, cmap=pylab.cm.binary, vmin=0, vmax=1)
	print "Writing '%s'" % (occupancy_fig,)
	pylab.savefig(occupancy_fig)
	pylab.show()

def plot_incidence(histories, incidence_fig):
	pylab.title("Occupied patches over time")
	pylab.ylabel("Occupied patches")
	pylab.xlabel("Time")
	for history in histories:
		pylab.plot(history / float(width*height))
	print "Writing '%s'" % (incidence_fig,)
	pylab.savefig(incidence_fig)
	pylab.show()

def plot():
	print "Plotting results for the simulations.."

	plot_habitat("habitat.pdf")

	for (index, s) in enumerate(spoms):
		plot_occupancy(s, "occupancy%i.pdf" % (index,))

	plot_incidence(histories, "incidence.pdf")

if do_plot:
	import pylab
	plot()
