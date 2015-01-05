# OVERVIEW
#
# A simulator script for working with the SPOM model.
# Input:
#  - habitat matrix
#  - list of species
#
# Optional input:
#  - region stencil file (label each habitat cell with a region ID)
#  - initial occupancy matrices for each species
# 
# Output:
# - occupancy matrices for each species 
# - regional occupancy statistics (how many patches in each region were occupied)
#
# USAGE
# 
# Run the script with parameter '-h' to get the list of options.
# To increase the level of logging run the script with parameters '-v' or '-vv' (for more detail).
#
# FILE FORMATS
#
# Habitat, stencil, and occupancy matrices are given as CSV files:
#  - comma ',' is used as a separator value
#  - row i, col j correspond to cell (i,j) in the matrix
# 
# Species are given as a CSV file. Each row corresponds to a species and the parameters
# are given in the following order:
#  - phenotype (phi)
#  - colonization rate (c)
#  - extinction rate (e)
#  - inverse of dispersal distance (alpha)
#  - niche width (gamma)

import scipy
import scipy.stats
import csv
import sys
import argparse
import os
import time
from multiprocessing import Pool

from collections import defaultdict

import model.grid
from model.grid import Species, Grid

VERBOSITY = 0

class InputData(object):
    def __init__(self, habitat, species_parameters, stencil, locations):
        self.habitat = habitat
        self.h, self.w = habitat.shape
        self.stencil = stencil
        self.species_parameters = species_parameters
        self.locations = locations
        self.species = [Species(*s) for s in self.species_parameters]

    def count_occupied_in_regions(self, occupancy):
        """ Count the number of occupied patches for each region specified in the stencil.
            If stencil is None, the whole grid is treated as a single region. """
        res = []

        if self.stencil is None:
            res.append(occupancy.sum())
        else:
            for r in scipy.unique(self.stencil):
                X = [self.stencil == r]*occupancy
                occurences = X.sum()
                assert occurences <= self.number_of_cells(r)
                res.append(occurences)

        return res

    def number_of_regions(self):
        if self.stencil is None:
            return 1
        else:
            return len(scipy.unique(self.stencil))

    def number_of_cells(self, region):
        if self.stencil == None:
            return self.h*self.w
        return (self.stencil.flatten() == region).sum()

    def create_grid(self, i):
        s = self.species[i]
        g = Grid(self.h, self.w, s)
        g.initialize()
        g.calculate_fitness_from_habitat(self.habitat)
        if len(args.stochasticity) > 0:
            mean, stddev, wf = args.stochasticity
            log("Enabling regional stochasticity with parameters: mean=%f, stddev=%f, weight factor=%f" % tuple(args.stochasticity), 1)
            g.enable_regional_stochasticity(mean, stddev, wf)

        # If a location matrix was specified, initialize the occupancy matrix with that,
        # otherwise populate everything.
        if i < len(self.locations):
            g.occupancy[:] = self.locations[i]
        else:
            g.occupancy[:] = 1
            
        return g

def read_matrix(csv_filename):
    """ Read a matrix from a CSV file"""

    log("Reading from '%s'" % (csv_filename,), 2)
    f = open(csv_filename)
    csv_reader = csv.reader(f)
    rows = []
    for row in csv_reader:
        rows.append(map(float, row))

    return scipy.array(rows)

def write_matrix(matrix, csv_filename):
    """ Write a matrix to a CSV file """
    log("Writing to '%s'" % csv_filename, 2)
    f = open(csv_filename, 'w')
    for row in matrix:
        f.write(",".join(map(str, row)))
        f.write("\n")
    f.close()

def read_input():
    """ Parse command line parameters """
    def assert_size(A,B, name1, name2):
        assert A.shape == B.shape, "%s and %s matrices have different shape: %s != %s" % (name1, name2, A, B)

    log("Reading input files..", 1)

    habitat = read_matrix(args.habitat)
    log("Landscape size %ix%i" % habitat.shape, 1)

    species = read_matrix(args.species)
    log("Loaded %i species" % (len(species,)), 1)

    stencil = None
    if args.stencil is not None:
        stencil = read_matrix(args.stencil) 
        assert_size(habitat, stencil, "habitat", "stencil")

    locations = [read_matrix(fname) for fname in args.locations]

    for i, l in enumerate(locations):
        assert_size(habitat, l, "habitat", "location %i" % (i,))
 
    return InputData(habitat, species, stencil, locations)

def initialize():
    """ Initialize the SPOM library """
    seed = int(time.time()*1000) if args.seed is None else args.seed
    log("Initializing the SPOM library (fft threads=%i, seed=%i)" % (args.fft_threads, seed), 1)
    model.grid.initialize_model(args.fft_threads, seed)


def sim_stencil_helper(p):
    i, data, time = p
    g = data.create_grid(i)
    hist = defaultdict(list)

    log("Simulating species %i" % i, 2)
    ps = g.occupancy.flatten().sum()
    for t in xrange(time+1):
        res = data.count_occupied_in_regions(g.occupancy)

        for region, patches in enumerate(res):
            hist[region].append(patches) 

        if ps == 0:
            log("Species %i EXTINCT at %i" % (i,t), 2)
            break # Extinct, no need to simulate

        ps = g.step()

    return hist, g.occupancy


def simulate_stencil(data, time):
    """ Compute number of occupied patches for each species in each region """
    # TODO This is a very slow and wasteful implementation
    species = len(data.species)
    regions = data.number_of_regions()
    total_hist = defaultdict(list)

    log("Using %i threads for simulation" % (args.threads,), 1)
    pool = Pool(args.threads)
    res_list = pool.map(sim_stencil_helper, [(i,data,time) for i in xrange(species)])

    assert len(res_list) == species, "Did not get data from all species"
    hists, occs = zip(*res_list)

    for i in xrange(len(hists)):
        for k in xrange(len(hists[i])):
            total_hist[(i,k)] = hists[i][k]

    log("Simulations complete")
    
    return total_hist, occs

def write_occupancy(occs):
    prefix = os.path.join(args.output_dir, args.prefix)
    # Occupancy matrices for each species 
    if args.output_locations:
        for i, om in enumerate(occs):
            write_matrix(om, prefix+'occupancy-%i.csv' % (i,))

def write_hist(hist, fname):
    log("Writing history to '%s'" % (fname,), 2)
    f = open(fname, 'w')
    for ((i,r), l) in sorted(hist.items()):
        line = [i,r] + l
        f.write(",".join(map(str,line)))
        f.write("\n")

def log(msg, lvl=1):
    timestamp = "%.3f" % (time.time() - START_TIME)
    if VERBOSITY >= lvl:
       print "%s -- %s" % (timestamp, msg)
       
def main():
    if args.output_dir is None:
        log("*** NO output files as output directory is not set ***", 0)

    initialize()
    data = read_input()

    r = data.number_of_regions()
    log("%i regions specified" % (r,), 1)
    for i in xrange(r):
        log("Region %i has %i cells" % (i, data.number_of_cells(i)), 2)
    
    for (i,s) in enumerate(data.species):
        log("Species %i: %s" % (i,s), 2)


    hist, occs = simulate_stencil(data, args.steps)

    if args.output_dir is not None:
        write_occupancy(occs)
        if hist is not None:
            fname = os.path.join(args.output_dir, args.prefix+"history.csv")
            write_hist(hist, fname)

    log("Finished", 0)

# --- argument parsing and calling of main ---

def make_parser():
    parser = argparse.ArgumentParser(description='Run a SPOM simulation from given input data')
    parser.add_argument('habitat', help='habitat matrix csv file')   
    parser.add_argument('species', help='species parameters csv file')
    parser.add_argument('--stencil', '-s', help='the region stencil csv file')
    parser.add_argument('--seed', type=int, help='the seed for the random number generator')
    parser.add_argument('--locations', '-l', default=[], nargs='+', help='initial occupancy matrices for species. if not specified, all species are present at each patch')
    parser.add_argument('--fft_threads', '-ft', type=int, default=1, help='number of threads allocated to FFTW')
    parser.add_argument('--threads', '-t', type=int, default=2, help='number of threads for parallel map')
    parser.add_argument('--steps', type=int, default=1, help='the number of simulation steps')
    parser.add_argument('--output_dir', '-od', help='write result files into the output directory')
    parser.add_argument('--output_locations', '-ol', action='store_true', default=False, help='output resulting occupancy matrices')
    parser.add_argument('--prefix', '-p', default="", help='prefix for output files')
    parser.add_argument('--stochasticity', type=float, nargs=3, default=[], help='enable regional stochasticity. Give mean, stddev, weight factor as parameters.')
    parser.add_argument('--verbose', '-v', action='count', default=0)
    return parser.parse_args()
 

START_TIME = time.time()
if __name__ == '__main__':
    args = make_parser()
    VERBOSITY = args.verbose
    main()
