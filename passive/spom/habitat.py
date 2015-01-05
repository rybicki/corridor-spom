"""
 Generate random landscapes with spatially auto-correlated habitat types.
"""

import numpy as np
import math

class HabitatGenerator(object):
    """ Generates spatially auto-correlated landscapes. The habitat matrix generated
        is a HxW matrix with values in [0,1]. 

        The generator is given the weight factor as parameter.

        After initialization the method generate() needs to be called to create the habitat
        matrix self.grid."""

    def __init__(self, height, width, weight_factor, random_shifts=True):
        self.height = int(height)
        self.width = int(width)
        assert self.height > 0 and self.width > 0, "Height and width must be positive"
        self.shape = (self.height, self.width)

        min_dim = min(self.shape)
        self.scale = int( math.log(min_dim, 2) ) + 1

        self.weight_factor = float(weight_factor)

        # Use randomized origo shifting?
        self.random_shifts = random_shifts

        # Give the zeroth level (i.e., the whole grid) zero weight, first level 1.0 
        self.__make_weights([0.0, 1.0])
        
    def generate(self):
        """ Generates a new habitat matrix and returns it. """
        self.grid = np.zeros( (self.shape) )

        # Generate origo shifts for each level before hand
        self.__make_shifts()

        # Generates autocorrelatet values in [0,1]
        self.__generate(0, (0,0), (self.shape))

        retval = self.grid
        self.grid = None # The object does not need to worry about the grid anymore
        return retval


    def __make_weights(self, initial_weights):
        self.weights = [] 
        self.weights.extend(initial_weights)
        
        for i in xrange(len(initial_weights), self.scale):
            new_weight = self.weights[i-1]/self.weight_factor
            self.weights.append(new_weight)

        # Normalize weights such that the sum of weights is 1.
        self.weights = np.array(self.weights)
        self.weights /= sum(self.weights)

    def __make_shifts(self):
        levels = len(self.weights)

        if not self.random_shifts:
            self.shifts = [(0,0)] * levels
            return

        # Randomize origo shifts
        self.shifts = []
        h,w = self.grid.shape

        for i in xrange(levels):
            dy = np.random.random_integers(-h, h)
            dx = np.random.random_integers(-w, w)
            self.shifts.append( (dy,dx) )

    def __generate(self, k, tl, br):
        y0, y1 = tl[0], br[0]
        x0, x1 = tl[1], br[1]

        vals = np.random.uniform() * self.weights[k]
        self.__update((y0,y1), (x0,x1), k, vals)

        # Recurse if maximum depth not exceeded yet
        if k+1 < self.scale:
            for (new_top_left, new_bottom_right) in self.__split_to_quads(tl, br):
                self.__generate(k+1, new_top_left, new_bottom_right)

    def __update(self, rangey, rangex, k, term):
        """ Update a region """
        minx, maxx = rangex
        miny, maxy = rangey
        modx, mody = self.grid.shape
        dy, dx = self.shifts[k]
        for i in xrange(miny, maxy):
            for j in xrange(minx, maxx):
                x = (j + dx) % modx
                y = (i + dy) % mody
                self.grid[y,x] += term

    def __split_to_quads(self, top_left, bottom_right):
        """ Split the square into four pieces. 
           The squares are returned as a list of
           top-left, bottom-right co-ordinate pairs."""
        x0,y0 = top_left
        x1,y1 = bottom_right
        # recursively split in to quads
        # a-b--
        # |A|B|
        # c-d-e
        # |C|D|
        # --f-g
        xx = x0 + (x1-x0)/2
        yy = y0 + (y1-y0)/2
        a = (x0, y0)
        b = (xx, y0)
        c = (x0, yy)
        d = (xx, yy)
        e = (x1, yy)
        f = (xx, y1)
        g = (x1, y1)

        return [ (a,d), (b,e), (c,f), (d,g) ]

def generate_habitat(height, width, weight_factor, random_shifts=True):
    """ Generates an autocorrelated habitat matrix of dimensions 'height' x 'width'. """
    gen = HabitatGenerator(height, width, weight_factor, random_shifts)
    return gen.generate()
