#@+leo-ver=4-thin
#@+node:gcross.20090813095726.2343:@thin sample.py
#@@language python
#@@tabwidth -4

import unittest
from paycheck import with_checker, positive_float, non_negative_float, irange, frange, unit_interval_float
from numpy import array, zeros, double, float64, isfinite, isfortran
from numpy.linalg import norm
from numpy.random import rand
from random import randint, random
import vpi.fortran as vpif

#@+others
#@+node:gcross.20090813095726.2344:sample_scheme1
class sample_scheme1(unittest.TestCase):
    #@    @+others
    #@+node:gcross.20090813095726.2345:test_move_type_within_bound
    @with_checker
    def test_move_type_within_bound(self,
            p1 = unit_interval_float, p2 = unit_interval_float, p3 = unit_interval_float,
            d1 = positive_float, d2 = positive_float, d3 = positive_float,
            n_slices = irange(2,5), n_particles = irange(2,5), n_dimensions = irange(2,5),
            lam = float,
        ):
        q0 = array(rand(n_slices,n_particles,n_dimensions),order='Fortran')
        q1 = array(rand(n_slices,n_particles,n_dimensions),order='Fortran')
        move_type_probabilities = array([p1,p2,p3])/(p1+p2+p3)
        move_type_differentials = array([d1,d2,d3])
        dM = randint(1,n_slices)
        swap_dimension_low = randint(1,n_dimensions)
        swap_dimension_high = randint(swap_dimension_low,n_dimensions)
        _,_,_,move_type = vpif.sample.sample_scheme1(q0,q1,move_type_probabilities,move_type_differentials,dM,lam,swap_dimension_low,swap_dimension_high)
        self.assert_(move_type >= 1)
        self.assert_(move_type <= 3)
    #@nonl
    #@-node:gcross.20090813095726.2345:test_move_type_within_bound
    #@+node:gcross.20090813095726.2347:test_move_type_fits_distribution
    @with_checker
    def test_move_type_fits_distribution(self,
            desired_move_type = irange(1,3),
            d1 = positive_float, d2 = positive_float, d3 = positive_float,
            n_slices = irange(2,5), n_particles = irange(2,5), n_dimensions = irange(2,5),
            lam = float,
        ):
        q0 = array(rand(n_slices,n_particles,n_dimensions),dtype='d',order='Fortran')
        q1 = array(rand(n_slices,n_particles,n_dimensions),dtype='d',order='Fortran')
        move_type_probabilities = zeros((3),dtype=double)
        move_type_probabilities[desired_move_type-1] = 1
        move_type_differentials = array([d1,d2,d3])
        dM = randint(1,n_slices)
        swap_dimension_low = randint(1,n_dimensions)
        swap_dimension_high = randint(swap_dimension_low,n_dimensions)
        _,_,_,move_type = vpif.sample.sample_scheme1(q0,q1,move_type_probabilities,move_type_differentials,dM,lam,swap_dimension_low,swap_dimension_high)
        self.assert_(move_type == desired_move_type)
    #@nonl
    #@-node:gcross.20090813095726.2347:test_move_type_fits_distribution
    #@-others
#@-node:gcross.20090813095726.2344:sample_scheme1
#@-others

tests = [
    sample_scheme1
    ]

if __name__ == "__main__":
    unittest.main()
#@nonl
#@-node:gcross.20090813095726.2343:@thin sample.py
#@-leo
