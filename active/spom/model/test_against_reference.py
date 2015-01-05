from nose.tools import *
import scipy
import cPickle
import pprint

from test_commons import *

import grid_model as gm 
from model import Species

def simulate(seed, height, width, species, habitat, steps):
    # Now create a normal simulation with same seed
    gm.init_fg(2, seed)
    g = gm.create_grid_from_struct(height, width, species)
    gm.initialize_simulation(g) # build the kernel
    gm.calculate_fitness(g, habitat)
    
    occ = gm.give_occupancy_matrix(g)
    occ[:] = 1
    gm.many_steps(g, steps)

    occ = scipy.copy( gm.give_occupancy_matrix(g) )
    fit = scipy.copy( gm.give_fitness_matrix(g) )
    gm.finalize_fg()
    return occ, fit

def dump_reference(seed, fname):
    print "Writing a reference file with seed '%i' to file '%s'" % (seed, fname)
    dim = scipy.random.randint(4, 128, 2)
    h,w = map(int, dim)
    habitat = scipy.random.random(h*w)
    habitat.resize( (h,w) )
    species = random_species(1, seed)[0]
    steps = scipy.random.randint(1, 20)

    ref_occ, ref_fit = simulate(seed, h, w, species.getModelStruct(), habitat, steps)
    ref_dict = {
        "height" : h,
        "width" : w,
        "occ" : ref_occ, 
        "fit" : ref_fit, 
        "species" : species.as_tuple(),
        "seed" : seed, 
        "habitat" : habitat, 
        "steps" : steps }

    pprint.pprint(ref_dict)

    cPickle.dump(ref_dict, open(fname, 'w'))

# The names of the reference files to read
reference_files = ["ref%i" % i for i in range(10)]
reference_files = ["test_refs/"+s for s in reference_files]

def make_references():
    for fname in reference_files:
        dump_reference(scipy.random.randint(1, 2**32), fname)

def check_reference(fname):
    print "Opening reference simulation '%s'" % fname
    ref = cPickle.load(open(fname, 'r')) 
    occ, fit = simulate(ref["seed"], ref["height"], ref["width"], Species(*ref["species"]).getModelStruct(), ref["habitat"], ref["steps"])
    assert occ.min() >= 0 and occ.max() <= 1, "Occupancy matrix needs to binary"
    compare_matrices(fit, ref["fit"], "Reference fitness mismatches. Why?")
    compare_matrices(occ, ref["occ"], "Reference occupancy mismatch!")

def test_against_references():
    """ Test against previously saved references to ensure that changes to the model do not break it"""
    for fname in reference_files:
        yield check_reference, fname

