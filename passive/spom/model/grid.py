"""
 This module wraps the 'grid_model' C library. A simulation of single species can be done with the Grid class.
"""

import numpy
import itertools as it
import grid_model as gm

__MODEL_INITIALIZED = False

def initialize_model(nthreads, seed):
    """ Call this before creating any Grid objects. """
    gm.init_fg(nthreads, seed)
    global __MODEL_INITIALIZED
    __MODEL_INITIALIZED = True
    # Automatically clean up
    import atexit
    atexit.register(finalize_model)

def is_model_initialized():
    global __MODEL_INITIALIZED
    return __MODEL_INITIALIZED

def finalize_model():
    if __MODEL_INITIALIZED:
        gm.finalize_fg()

class Species(object):
    """ Store the parameters associated to a single species"""

    def __init__(self, phenotype, colonization_rate, extinction_rate, alpha, bell_width):
        """ Create a new Species object. """
        self._phenotype = phenotype
        self._colonization_rate = colonization_rate
        self._extinction_rate = extinction_rate
        self._alpha = alpha
        self._fitness_bell_width = bell_width

    def getModelStruct(self):
        """ Return a grid_model.c module representation of species parameters. 
            This is the t_species_parameters struct defined in grid_model.h """
        sp = gm.create_species_parameters( self.phenotype, self.colonization_rate, self.extinction_rate, self.alpha, self.fitness_bell_width)
        return sp
    
    def as_tuple(self):
        """ This is for compatibility reasons. """
        return (self.phenotype, self.colonization_rate, self.extinction_rate, self.alpha, self.fitness_bell_width)

    def __getitem__(self, i):
        return self.as_tuple()[i]

    def __repr__(self):
        return "Species(phenotype=%f, colonization_rate=%f, extinction_rate=%f, alpha=%f, fitness_bell_width=%f)" % self.as_tuple()

    def __str__(self):
        return self.__repr__()

    def __eq__(self, other):
        return isinstance(other, self.__class__) and self.__dict__ == other.__dict__

    def __ne__(self, other):
        return not self.__eq__(other)

    phenotype = property(fget=lambda self: self._phenotype)
    colonization_rate = property(fget=lambda self: self._colonization_rate)
    extinction_rate = property(fget=lambda self: self._extinction_rate)
    alpha = property(fget=lambda self: self._alpha)
    fitness_bell_width = property(fget=lambda self: self._fitness_bell_width)


class NotInitializedException(Exception):
    pass

class Grid(object):
    """ The grid of all occupancy patches for a single species. """

    def __init__(self, height, width, species_parameters):
        """ Create an uninitialized Grid. Be sure to call initialize_model() before hand. """
        self._species_parameters = species_parameters
        self._t_grid = gm.create_grid_from_struct(height, width, self.species_parameters.getModelStruct())
        self._height = height
        self._width = width
        self._initialized = False

    def __del__(self):
        if self._t_grid != None:
            gm.free_grid(self._t_grid)

    def free(self):
        """ Free the resources used by this grid."""
        gm.free_grid(self._t_grid)
        self._t_grid = None
        self._initialized = False

    def set_occupancy(self, occ):
        """ Copy the contents of 'occ' to the grid's occupancy array. The array 'occ' should be of type that can be safely casted into numpy.uint8"""
        if occ.shape != (self.height, self.width):
            raise ValueError("Occupancy matrix has wrong dimensions")
        gm.copy_occupancy_from(self._t_grid, occ)

    def get_occupancy(self):
        return gm.give_occupancy_matrix(self._t_grid)

    def enable_regional_stochasticity(self, mean, stddev, weight_factor):
        """ Enable regional stochasticity with the given mean, standard deviation (sigma) and weight factor """
        p = gm.t_stochasticity_parameters()
        p.mean = mean
        p.stddev = stddev
        p.weight_factor = weight_factor
        gm.enable_regional_stochasticity(self._t_grid, p)

    def disable_regional_stochasticity(self):
        """ Disable regional stochasticity. """
        gm.disable_regional_stochasticity(self._t_grid)

    occupancy = property(fset=set_occupancy, fget=get_occupancy)
    fitness = property(fget=lambda self: gm.give_fitness_matrix(self._t_grid)) 
    species_parameters = property(fget=lambda self: self._species_parameters)
    height = property(fget=lambda self: self._height)
    width = property(fget=lambda self: self._width)

    def initialize(self):
        """ Initialize the grid (build the dispersal kernel etc.) """
        gm.initialize_simulation(self._t_grid);
        self._initialized = True

    def is_initialized(self):
        return is_model_initialized() and self._initialized

    def step(self):
        """ Perform a single step of the simulation. """
        if not self.is_initialized(): 
            raise NotInitializedException("Grid not initialized")
        return gm.simulate_step(self._t_grid)

    def many_steps(self, n):
        """ Perform 'n' steps of the simulation. """
        if not self.is_initialized(): 
            raise NotInitializedException("Grid not initialized")
        return gm.many_steps(self._t_grid, n)

    def calculate_fitness_from_habitat(self, habitat):
        """ Calculate the species fitness from the given habitat matrix. """
        if habitat.shape != (self.height, self.width):
            raise ValueError("Habitat size must match grid shape")
        gm.calculate_fitness(self._t_grid, habitat)

