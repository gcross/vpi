#@+leo-ver=4-thin
#@+node:gcross.20090806151612.1843:@thin test.py
#@@language python
#@@tabwidth -4

import unittest
from paycheck import with_checker, generator
import numpy
from numpy import array, zeros, all
from numpy.random import rand

class ArrayGenerator(generator.CollectionGenerator):
    def __init__(self, dimensions, dtype=numpy.double, num_calls=generator.NUM_CALLS):
        self.dimensions = dimensions
        self.dtype = dtype

    def next_value(self):
        return numpy.array(rand(*[random.randint(0, generator.LIST_LEN) for dummy in self.dimensions]),dtype=self.dtype)

class update_xij(unittest.TestCase):
    @with_checker(ArrayGenerator(2))
    def values_are_nonnegative(self,x):
        n_slices, n_particles = x.shape
        xij2 = zeros((n_slices,n_particles,n_particles))
        vpi.xij.update_xij(xij2,x,1,n_slices)
        self.assert_(all(xij2>=0))

if __name__ == "__main__":
    unittest.main()
#@-node:gcross.20090806151612.1843:@thin test.py
#@-leo
