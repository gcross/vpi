#@+leo-ver=4-thin
#@+node:gcross.20090806151612.1843:@thin xij.py
#@@language python
#@@tabwidth -4

import unittest
from paycheck import with_checker, positive_float, non_negative_float, irange
from numpy import zeros, double, array
from numpy.random import rand
from numpy.linalg import norm
import vpi.fortran as vpif

#@+others
#@+node:gcross.20090807144330.1678:wrap_around
class wrap_around(unittest.TestCase):

    @with_checker
    def test_values_within_correct_range(self,x = positive_float,period_length = positive_float(minimum_magnitude=1)):
        period_length = abs(period_length)
        result = vpif.xij.wrap_around(x,period_length)
        self.assert_(result <= period_length/2)
        self.assert_(result >= -period_length/2)
#@nonl
#@-node:gcross.20090807144330.1678:wrap_around
#@+node:gcross.20090807144330.1672:compute_xij
class compute_xij(unittest.TestCase):
    #@    @+others
    #@+node:gcross.20100922111545.1737:test_values_are_nonnegative
    @with_checker
    def test_values_are_nonnegative(self,n_slices = irange(1,5), n_particles = irange(2,5), n_dimensions = irange(1,4)):
        x = array(rand(n_dimensions,n_particles,n_slices),dtype=double,order='Fortran')
        xij2 = vpif.xij.compute_xij(x)
        self.assertEqual(xij2.shape,(n_particles,n_particles,n_slices))
        self.assert_((xij2>=0).all())
    #@nonl
    #@-node:gcross.20100922111545.1737:test_values_are_nonnegative
    #@+node:gcross.20100922111545.1738:test_values_are_correct
    @with_checker
    def test_values_are_correct(self,n_slices = irange(1,5), n_particles = irange(2,5), n_dimensions = irange(1,4), s = int, i = int, j = int):
        x = array(rand(n_dimensions,n_particles,n_slices),dtype=double,order='Fortran')
        xij2 = vpif.xij.compute_xij(x)
        self.assertEqual(xij2.shape,(n_particles,n_particles,n_slices))
        s = s % n_slices
        i = i % n_particles
        j = j % n_particles
        self.assertAlmostEqual(xij2[j,i,s],norm(x[:,i,s]-x[:,j,s])**2)
    #@nonl
    #@-node:gcross.20100922111545.1738:test_values_are_correct
    #@-others
#@-node:gcross.20090807144330.1672:compute_xij
#@+node:gcross.20090807144330.1680:compute_xij_pbc
class compute_xij_pbc(unittest.TestCase):
    #@    @+others
    #@+node:gcross.20090818081913.1350:test_values_are_nonnegative
    @with_checker
    def test_values_are_nonnegative(self, n_slices = irange(1,5), n_particles = irange(1,2), n_dimensions = irange(1,4), period_length = positive_float(1)):
        x = array(rand(n_dimensions,n_particles,n_slices) % (period_length/2),dtype=double,order='Fortran')
        xij2 = vpif.xij.compute_xij_pbc(x,period_length)
        self.assertEqual(xij2.shape,(n_particles,n_particles,n_slices))
        self.assert_((xij2>=0).all())
    #@nonl
    #@-node:gcross.20090818081913.1350:test_values_are_nonnegative
    #@+node:gcross.20090818081913.1351:test_values_are_within_range
    @with_checker
    def test_values_are_within_range(self, n_slices = irange(1,5), n_particles = irange(1,2), n_dimensions = irange(1,4), period_length = positive_float(1)):
        x = array(rand(n_dimensions,n_particles,n_slices) % (period_length/2),dtype=double,order='Fortran')
        xij2 = vpif.xij.compute_xij_pbc(x,period_length)
        self.assertEqual(xij2.shape,(n_particles,n_particles,n_slices))
        self.assert_((xij2<=period_length/2).all())
    #@nonl
    #@-node:gcross.20090818081913.1351:test_values_are_within_range
    #@+node:gcross.20090818081913.1352:test_values_are_correct
    @with_checker
    def test_values_are_correct(self, n_slices = irange(1,5), n_particles = irange(1,2), n_dimensions = irange(1,4), s = int, i = int, j = int, period_length = positive_float(1)):
        x = array(rand(n_dimensions,n_particles,n_slices) % (period_length/2),dtype=double,order='Fortran')
        xij2 = vpif.xij.compute_xij_pbc(x,period_length)
        self.assertEqual(xij2.shape,(n_particles,n_particles,n_slices))
        s = s % n_slices
        i = i % n_particles
        j = j % n_particles
        self.assertAlmostEqual(xij2[j,i,s],(norm(x[:,i,s]-x[:,j,s])%(period_length/2))**2)
    #@nonl
    #@-node:gcross.20090818081913.1352:test_values_are_correct
    #@-others

#@-node:gcross.20090807144330.1680:compute_xij_pbc
#@-others

tests = [
    wrap_around,
    compute_xij,
    compute_xij_pbc,
    ]

if __name__ == "__main__":
    unittest.main()
#@nonl
#@-node:gcross.20090806151612.1843:@thin xij.py
#@-leo
