#@+leo-ver=4-thin
#@+node:gcross.20090807144330.2151:@thin thermalize.py
#@@language python
#@@tabwidth -4

import unittest
from paycheck import with_checker
from paycheck.generator import positive_float, non_negative_float, irange, frange
from numpy import array, zeros, double, float64, isfinite, int32
from numpy.linalg import norm
from numpy.random import rand
from random import randint, random
import vpi.fortran as vpif

#@+others
#@+node:gcross.20090807144330.2152:accept_path
class accept_path(unittest.TestCase):
    @with_checker
    def test_always_accepts_minimizing_move(self, p1 = positive_float, p2 = positive_float):
        self.assert_(
            (vpif.thermalize.accept_path(p1,p2) == 1) or
            (p1 > p2)
        )
#@nonl
#@-node:gcross.20090807144330.2152:accept_path
#@+node:gcross.20090813095726.2348:thermalize_path
class thermalize_path(unittest.TestCase):
    #@    @+others
    #@+node:gcross.20090813095726.2349:test_rejection_case
    @with_checker(number_of_calls=10)
    def test_rejection_case(self,
            p1 = non_negative_float, p2 = non_negative_float, p3 = non_negative_float,
            d1 = float, d2 = float, d3 = float,
            n_trials = irange(1,10),
            center_slice = irange(3,11,2),
            n_particles = irange(1,10),
            lam = float,
        ):
        N_rotating_particles = randint(1,n_particles)
        n_slices = center_slice * 2
        n_dimensions = 3
        x = zeros((n_dimensions,n_particles,n_slices),dtype=double,order='Fortran')
        xij2 = zeros((n_particles,n_particles,n_slices),dtype=double,order='Fortran')
        U = zeros((n_particles,n_slices),dtype=double,order='Fortran')
        gradU2 = zeros((n_slices,),dtype=double,order='Fortran')
        particle_number = randint(1,n_particles)
        move_type_probabilities = array([p1,p2,p3])/(p1+p2+p3)
        move_type_differentials = array([d1,d2,d3])
        dM = randint(1,n_slices)
        swap_dimension_low = randint(1,n_dimensions)
        swap_dimension_high = randint(swap_dimension_low,n_dimensions)
        def null_func(*args):
            return double(0.0)
        def compute_potential(x,xij2,n_slices,n_particles,n_dimensions):
            if((abs(x)==0).all()):
                return zeros((n_particles,n_slices)), zeros((n_slices,)), False
            else:
                return zeros((n_particles,n_slices)), zeros((n_slices,)), True
        def greens_function(*args):
            return double(0.0)
        slice_move_attempted_counts, slice_move_accepted_counts = [zeros((n_slices,),dtype='i',order='Fortran') for dummy in xrange(2)]
        move_type_attempted_counts, move_type_accepted_counts = [zeros((3,),dtype='i',order='Fortran') for dummy in xrange(2)]
        vpif.thermalize.thermalize_path(
            x,xij2,
            U,gradU2,
            n_trials,
            move_type_probabilities,move_type_differentials,
            dM,lam,swap_dimension_low,swap_dimension_high,
            slice_move_attempted_counts,move_type_attempted_counts,slice_move_accepted_counts,move_type_accepted_counts,
            compute_potential, null_func,greens_function,
        )
        self.assert_((x==0.0).all())
        self.assert_((xij2==0.0).all())
    #@nonl
    #@-node:gcross.20090813095726.2349:test_rejection_case
    #@-others
#@-node:gcross.20090813095726.2348:thermalize_path
#@-others

tests = [
    accept_path,
    thermalize_path
    ]

if __name__ == "__main__":
    unittest.main()
#@nonl
#@-node:gcross.20090807144330.2151:@thin thermalize.py
#@-leo
