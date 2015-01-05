from nose.tools import *
import time
import scipy

import model
from model import Species, Grid

class TestSpecies:
    """ Unit test for model.Species class """
    def setUp(self):
        self.species_inst = Species(0.1, 0.2, 0.3, 0.4, 0.5)

    def tearDown(self):
        self.species_inst = None

    def test_eq(self):
        inst1 = Species(0.1,0.2,0.3,0.4,0.5)
        inst2 = Species(0.1,0.2,0.3,0.4,0.5)
        assert inst1 is not inst2
        assert inst1 == inst2
        assert (inst1 == inst2) == (not (inst1 != inst2)), "__eq__ and __ne__ should be complementary"
        assert not (inst1 != inst2), "Contests are equal so not neq should be true"
        inst3 = Species(0.9, 0.2, 0.3, 0.4, 0.5)
        assert not (inst1 == inst3)
        assert inst1 != inst3
        
    def test_tuple(self):
        t = (0.1, 0.2, 0.3, 0.4, 0.5)
        s = Species(*t)
        for i in xrange(len(t)):
            assert s[i] == t[i]

    @raises(AttributeError)
    def test_immutable_phenotype(self):
       self.species_inst.phenotype = 0.5

    @raises(AttributeError)
    def test_immutable_alpha(self):
       self.species_inst.alpha = 0.5

    @raises(AttributeError)
    def test_immutable_extinction(self):
       self.species_inst.extinction_rate = 0.5

    @raises(AttributeError)
    def test_immutable_colonization(self):
       self.species_inst.colonization_rate = 0.5

    @raises(AttributeError)
    def test_immutable_fitness(self):
       self.species_inst.fitness_bell_width = 0.5

MTHREADS = 2

class TestGrid:
    """ Unit test for model.Grid class """

    def setUp(self):
        model.initialize_model(MTHREADS, int(time.time()))

    def tearDown(self):
        model.finalize_model()

    @nottest
    def generate_test_input(self):
        shapes = [ (1,1), (2,2), (32, 32), (7, 13), (17, 32) ]
        species = [ Species(0.1, 0.2, 0.3, 0.4, 0.5), Species(0.5, 0.5, 0.5, 0.5, 0.5) ]
        for sh in shapes:
            for sp in species:
                yield sh, sp

    def check_create(self, shape, species):
        g = Grid(shape[0], shape[1], species)
        assert g.height == shape[0]
        assert g.width == shape[1]
        assert g.species_parameters == species
        assert not g.is_initialized()
        g.initialize()
        assert g.is_initialized()

    def test_create(self):
        for (shape, species) in self.generate_test_input():
            yield self.check_create, shape, species

    def check_free(self, shape, species):
        g = Grid(shape[0], shape[1], species)
        g.free()
        assert not g.is_initialized()
        del g
        
        g = Grid(shape[0], shape[1], species)
        g.initialize()
        g.free()
        assert not g.is_initialized()
        del g

    def test_free(self):
        for (shape, species) in self.generate_test_input():
            yield self.check_create, shape, species

    def check_occupancy_accessors(self, shape, species):
        g = Grid(shape[0], shape[1], species)
        occ = scipy.ones( (shape[0], shape[1]), dtype=scipy.uint8)
        occ[0,0] = 0
        g.occupancy = occ
        assert g.occupancy is not occ, "The grid should not just copy a reference but the contents"
        assert (g.occupancy == occ).all(), "The grid should copy the occupancy contents correctly"
        occ[0,0] = 1
        assert not (g.occupancy == occ).all(), "Changes to occ should not appear in the grid"

    def test_occupancy_accessors(self):
        for (shape, species) in self.generate_test_input():
            yield self.check_occupancy_accessors, shape, species

    @raises(model.NotInitializedException)
    def check_step_not_initialized(self, shape, species):
        g = Grid(shape[0], shape[1], species)
        g.step()

    @raises(model.NotInitializedException)
    def check_many_steps_not_initialized(self, shape, species):
        (shape, species) = self.generate_test_input().next()
        g = Grid(shape[0], shape[1], species)
        g.many_steps(10)

    def test_steps_not_initialized(self):
        for (shape, species) in self.generate_test_input():
            yield self.check_step_not_initialized, shape, species
            yield self.check_many_steps_not_initialized, shape, species

    def check_step(self, g):
        ret = g.step()
        assert ret >= 0, "step should return the number of alive patches"

    def check_many_steps(self,g):
        for i in xrange(10):
            ret = g.many_steps(i)
            assert len(ret) == i, "There should be patch count for each time step"
            assert all([r >= 0 for r in ret]), "There cannot be negative number of alive patches"

    def test_extinct_stay_extinct(self):
        for (shape, species) in self.generate_test_input():
            g = Grid(shape[0], shape[1], species)
            habitat = scipy.ones( (shape) ) * 0.5
            g.calculate_fitness_from_habitat(habitat)
            g.occupancy[:] = 0
            g.initialize()
            assert all( [i == 0 for i in g.many_steps(25)] ), "Extinct species reappear!"
            g.enable_regional_stochasticity(0, 0.8, 1.0)
            assert all( [i == 0 for i in g.many_steps(25)] ), "Extinct species reappear!"



    def test_simulation_steps(self):
        for (shape, species) in self.generate_test_input():
            g = Grid(shape[0], shape[1], species)
            habitat = scipy.ones( (shape) ) * 0.5
            g.calculate_fitness_from_habitat(habitat)
            g.occupancy[:] = 1
            g.initialize()
            yield self.check_step, g
            yield self.check_many_steps, g

    def test_calculate_fitness(self):
        for (shape, species) in self.generate_test_input():
            g = Grid(shape[0], shape[1], species)
            habitat = scipy.ones( (shape) ) * 0.5
            g.calculate_fitness_from_habitat(habitat)
            habitat = scipy.ones( (shape) ) * species.phenotype
            g.calculate_fitness_from_habitat(habitat)
            assert (g.fitness == scipy.ones( (shape) )).all(), "Optimal fitness should be 1"

    def test_enable_regional_stochasticity(self):
        for (shape, species) in self.generate_test_input():
            g = Grid(shape[0], shape[1], species)
            habitat = scipy.ones( (shape) ) * 0.5
            g.calculate_fitness_from_habitat(habitat)
            g.occupancy[:] = 1
            g.initialize()

            assert g._t_grid.use_regional_stochasticity == 0, "use_regional_stochasticity flag should be initially false"
            assert g._t_grid.stochasticity_grid is None, "stochasticity_grid should be disabled initially"

            g.enable_regional_stochasticity(0, 0.8, 1.0)
            assert g._t_grid.use_regional_stochasticity > 0, "t_grid regional_stochasticity flag not set"
            assert g._t_grid.stochasticity_grid is not None
            yield self.check_many_steps, g

            g.disable_regional_stochasticity()

            assert g._t_grid.use_regional_stochasticity == 0, "use_regional_stochasticity flag should be false"
            assert g._t_grid.stochasticity_grid is None, "stochasticity_grid should be disabled"

