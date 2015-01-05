"""
 Run a single species simulation and plot each time step separately.
"""
import time
import numpy
import pylab

import model.grid
from model.grid import Species, Grid
import habitat

# The size of the simulation grid
height, width = 128, 128

# Simulation time
time_steps = 10

# Random number generator seed for the simulation
seed = int(time.time())

# Create a random autocorrelated habitat matrix. 
habitat_matrix = habitat.generate_habitat(height, width, 1.5, False)

model.grid.initialize_model(1, seed)

# Initialize a SPOM grid
species = Species(0.8, 1.0, 0.2, 0.5, 0.1) 
g = Grid(height, width, species)
g.initialize()
g.calculate_fitness_from_habitat(habitat_matrix) 
g.occupancy[:] = 1 # Initialize the occupancy matrix to contain all 1s

# Run the simulation and plot the occupancy matrix at each time step
for i in xrange(time_steps+1):
	pylab.title("Step %i/%i" % (i,time_steps))
	pylab.imshow(g.occupancy, cmap=pylab.cm.binary, vmin=0, vmax=1)
	pylab.show()
	g.step()
