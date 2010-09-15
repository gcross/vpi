#@+leo-ver=4-thin
#@+node:gcross.20090819093822.1403:@thin lattice.py
import unittest
from paycheck import *
from numpy import *
from numpy.random import rand
import __builtin__
import vpi.fortran as vpif

#@+others
#@+node:gcross.20090819093822.1404:class make_lattice
class make_lattice(unittest.TestCase):
    #@    @+others
    #@+node:gcross.20090819093822.1405:test_within_size
    @with_checker
    def test_within_size(self,
            n_slices=irange(1,5),
            n_particles=irange(1,10),
            n_dimensions=irange(1,4),
            size=unit_interval_float,
            ):
        lattice = vpif.lattice.make_lattice(size,n_slices,n_particles,n_dimensions)
        self.assert_((lattice>=-size/2).all())
        self.assert_((lattice<=size/2).all())
    #@nonl
    #@-node:gcross.20090819093822.1405:test_within_size
    #@+node:gcross.20090819093822.1409:test_finite
    @with_checker
    def test_finite(self,
            n_slices=irange(1,5),
            n_particles=irange(1,10),
            n_dimensions=irange(1,4),
            size=float,
            ):
        lattice = vpif.lattice.make_lattice(size,n_slices,n_particles,n_dimensions)
        self.assert_(isfinite(lattice).all())
    #@nonl
    #@-node:gcross.20090819093822.1409:test_finite
    #@-others
#@-node:gcross.20090819093822.1404:class make_lattice
#@-others

tests = [
    make_lattice,
    ]

if __name__ == "__main__":
    unittest.main()
#@nonl
#@-node:gcross.20090819093822.1403:@thin lattice.py
#@-leo
