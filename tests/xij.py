#@+leo-ver=4-thin
#@+node:gcross.20090806151612.1843:@thin xij.py
#@@language python
#@@tabwidth -4

import unittest
from paycheck import with_checker, positive_float, non_negative_float, irange
from numpy import zeros, double, array
from numpy.random import rand
from numpy.linalg import norm
import vpi

#@+others
#@+node:gcross.20090807144330.1678:wrap_around
class wrap_around(unittest.TestCase):

    @with_checker
    def test_values_within_correct_range(self,x = positive_float,period_length = positive_float(minimum_magnitude=1)):
        period_length = abs(period_length)
        result = vpi.xij.wrap_around(x,period_length)
        self.assert_(result <= period_length/2)
        self.assert_(result >= -period_length/2)
#@-node:gcross.20090807144330.1678:wrap_around
#@+node:gcross.20090807144330.1672:update_xij
class update_xij(unittest.TestCase):

    @with_checker
    def test_values_are_nonnegative(self,n_slices = irange(1,5), n_particles = irange(2,5)):
        x = array(rand(n_slices,n_particles),dtype=double,order='Fortran')
        xij2 = zeros((n_slices,n_particles,n_particles),dtype=double,order='Fortran')
        vpi.xij.update_xij(xij2,x)
        self.assert_((xij2>=0).all())

    @with_checker
    def test_values_are_correct(self,n_slices = irange(1,5), n_particles = irange(2,5), s = int, i = int, j = int):
        x = array(rand(n_slices,n_particles),dtype=double,order='Fortran')
        xij2 = zeros((n_slices,n_particles,n_particles),dtype=double,order='Fortran')
        vpi.xij.update_xij(xij2,x)
        s = s % n_slices
        i = i % n_particles
        j = j % n_particles
        self.assertAlmostEqual(xij2[s,i,j],norm(x[s,i]-x[s,j])**2)
#@-node:gcross.20090807144330.1672:update_xij
#@+node:gcross.20090807144330.1680:update_xij_pbc
class update_xij_pbc(unittest.TestCase):
    #@    @+others
    #@+node:gcross.20090818081913.1350:test_values_are_nonnegative
    @with_checker
    def test_values_are_nonnegative(self, n_slices = irange(1,5), n_particles = irange(1,2), period_length = positive_float(1)):
        x = array(rand(n_slices,n_particles) % (period_length/2),dtype=double,order='Fortran')
        xij2 = zeros((n_slices,n_particles,n_particles),dtype=double,order='Fortran')
        vpi.xij.update_xij_pbc(xij2,x,period_length)
        self.assert_((xij2>=0).all())
    #@-node:gcross.20090818081913.1350:test_values_are_nonnegative
    #@+node:gcross.20090818081913.1351:test_values_are_within_range
    @with_checker
    def test_values_are_within_range(self, n_slices = irange(1,5), n_particles = irange(1,2), period_length = positive_float(1)):
        x = array(rand(n_slices,n_particles) % (period_length/2),dtype=double,order='Fortran')
        xij2 = zeros((n_slices,n_particles,n_particles),dtype=double,order='Fortran')
        vpi.xij.update_xij_pbc(xij2,x,period_length)
        self.assert_((xij2<=period_length/2).all())
    #@nonl
    #@-node:gcross.20090818081913.1351:test_values_are_within_range
    #@+node:gcross.20090818081913.1352:test_values_are_correct
    @with_checker
    def test_values_are_correct(self, n_slices = irange(1,5), n_particles = irange(1,2), s = int, i = int, j = int, period_length = positive_float(1)):
        x = array(rand(n_slices,n_particles) % (period_length/2),dtype=double,order='Fortran')
        xij2 = zeros((n_slices,n_particles,n_particles),dtype=double,order='Fortran')
        vpi.xij.update_xij_pbc(xij2,x,period_length)
        s = s % n_slices
        i = i % n_particles
        j = j % n_particles
        self.assertAlmostEqual(xij2[s,i,j],(norm(x[s,i]-x[s,j])%(period_length/2))**2)
    #@nonl
    #@-node:gcross.20090818081913.1352:test_values_are_correct
    #@-others

#@-node:gcross.20090807144330.1680:update_xij_pbc
#@+node:gcross.20090828201103.2129:hard_sphere_violation
class hard_sphere_violation(unittest.TestCase):

    def test_violated(self):
        xij2 = array([[[0,3],[3,0]],[[0,1],[1,0]]],dtype=double,order='Fortran')
        hard_sphere_radius_squared = 2
        self.assert_(vpi.xij.hard_sphere_violation(xij2,hard_sphere_radius_squared))

    def test_not_violated(self):
        xij2 = array([[[0,3],[3,0]],[[0,2.5],[2.5,0]]],dtype=double,order='Fortran')
        hard_sphere_radius_squared = 2
        self.assert_(not vpi.xij.hard_sphere_violation(xij2,hard_sphere_radius_squared))
#@-node:gcross.20090828201103.2129:hard_sphere_violation
#@-others

tests = [
    wrap_around,
    update_xij,
    update_xij_pbc,
    hard_sphere_violation,
    ]

if __name__ == "__main__":
    unittest.main()
#@-node:gcross.20090806151612.1843:@thin xij.py
#@-leo
