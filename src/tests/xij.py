#@+leo-ver=4-thin
#@+node:gcross.20090806151612.1843:@thin xij.py
#@@language python
#@@tabwidth -4

import unittest
from paycheck import with_checker
from numpy import zeros, double
from numpy.linalg import norm
from tests import particle_paths_type
import vpi

#@+others
#@+node:gcross.20090807144330.1678:wrap_around
class wrap_around(unittest.TestCase):

    @with_checker(float,float)
    def test_values_within_correct_range(self,x,period_length):
        period_length = abs(period_length)
        result = vpi.xij.wrap_around(x,period_length)
        self.assert_(result <= period_length/2)
        self.assert_(result >= -period_length/2)
#@-node:gcross.20090807144330.1678:wrap_around
#@+node:gcross.20090807144330.1672:update_xij
class update_xij(unittest.TestCase):

    @with_checker(particle_paths_type)
    def test_values_are_nonnegative(self,x):
        n_slices, n_particles = x.shape
        xij2 = zeros((n_slices,n_particles,n_particles),dtype=double,order='Fortran')
        vpi.xij.update_xij(xij2,x,1,n_slices)
        self.assert_((xij2>=0).all())

    @with_checker(particle_paths_type,int,int,int)
    def test_values_are_correct(self,x,s,i,j):
        n_slices, n_particles = x.shape
        xij2 = zeros((n_slices,n_particles,n_particles),dtype=double,order='Fortran')
        vpi.xij.update_xij(xij2,x,1,n_slices)
        s = s % n_slices
        i = i % n_particles
        j = j % n_particles
        self.assertAlmostEqual(xij2[s,i,j],norm(x[s,i]-x[s,j])**2)
#@-node:gcross.20090807144330.1672:update_xij
#@-others

tests = [
    wrap_around,
    update_xij,
    ]

if __name__ == "__main__":
    unittest.main()
#@-node:gcross.20090806151612.1843:@thin xij.py
#@-leo
