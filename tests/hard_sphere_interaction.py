#@+leo-ver=4-thin
#@+node:gcross.20090911091023.2483:@thin hard_sphere_interaction.py
import unittest
from paycheck import *
from numpy import *
from numpy.random import rand
import __builtin__
import vpi.fortran as vpif

#@+others
#@+node:gcross.20090828201103.2129:has_collision
class has_collision(unittest.TestCase):

    def test_violated(self):
        xij2 = array([[[0,3],[3,0]],[[0,1],[1,0]]],dtype=double,order='Fortran')
        hard_sphere_radius_squared = 2
        self.assert_(vpif.hard_sphere_interaction.has_collision(xij2,hard_sphere_radius_squared))

    def test_not_violated(self):
        xij2 = array([[[0,3],[3,0]],[[0,2.5],[2.5,0]]],dtype=double,order='Fortran')
        hard_sphere_radius_squared = 2
        self.assert_(not vpif.hard_sphere_interaction.has_collision(xij2,hard_sphere_radius_squared))
#@nonl
#@-node:gcross.20090828201103.2129:has_collision
#@+node:gcross.20090828095451.1461:compute_greens_function
class compute_greens_function(unittest.TestCase):

    @with_checker
    def test_finite(self,
            n_slices = irange(1,10),
            n_particles = irange(1,5),
            dt = float,
            hard_sphere_radius = float,
        ):
        xij2 = rand(n_slices,n_particles,n_particles)
        slice_start = randint(1,n_slices)
        slice_end = randint(slice_start,n_slices)
        particle_number = randint(1,n_particles)
        vpif.hard_sphere_interaction.compute_greens_function(
            xij2,
            dt,hard_sphere_radius,
            slice_start, slice_end,
            particle_number
        )
#@nonl
#@-node:gcross.20090828095451.1461:compute_greens_function
#@+node:gcross.20090828095451.1463:compute_greens_function2
class compute_greens_function2(unittest.TestCase):

    @with_checker
    def test_finite(self,
            n_slices = irange(1,10),
            n_particles = irange(1,5),
            n_dimensions = irange(1,5),
            lambda_ = float,
            dt = float,
            hard_sphere_radius = float,
        ):
        q = rand(n_slices,n_particles,n_dimensions)
        slice_start = randint(1,n_slices)
        slice_end = randint(slice_start,n_slices)
        particle_number = randint(1,n_particles)
        vpif.hard_sphere_interaction.compute_greens_function2(
            q,
            lambda_,dt,hard_sphere_radius,
            slice_start, slice_end,
            particle_number
        )
#@nonl
#@-node:gcross.20090828095451.1463:compute_greens_function2
#@-others

tests = [
    has_collision,
    compute_greens_function,
    compute_greens_function2,
    ]

if __name__ == "__main__":
    unittest.main()
#@nonl
#@-node:gcross.20090911091023.2483:@thin hard_sphere_interaction.py
#@-leo
