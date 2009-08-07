#@+leo-ver=4-thin
#@+node:gcross.20090806151612.1843:@thin test.py
#@@language python
#@@tabwidth -4

import unittest
from paycheck import with_checker, generator
import numpy
from numpy import array, zeros, all, double
from numpy.linalg import norm
from numpy.random import rand
import random
import vpi

# ------------------------------------------------------------------------------
# NumPy Support
# ------------------------------------------------------------------------------

class ArrayGenerator(generator.PayCheckGenerator):
    def __init__(self, example):
        self.shape = example.shape
        self.dtype = example.dtype

    def next(self):
        return numpy.array(rand(*[random.randint(1, generator.LIST_LEN) for _ in self.shape]),dtype=self.dtype)

generator.container_generators[type(zeros(()))] = ArrayGenerator

particle_paths_type = zeros((0,0),dtype=double)

class update_xij(unittest.TestCase):

    @with_checker(particle_paths_type)
    def test_xij_values_are_nonnegative(self,x):
        n_slices, n_particles = x.shape
        xij2 = zeros((n_slices,n_particles,n_particles),dtype=double,order='Fortran')
        vpi.xij.update_xij(xij2,x,1,n_slices)
        self.assert_(all(xij2>=0))

    @with_checker(particle_paths_type,int,int,int)
    def test_xij_values_are_correct(self,x,s,i,j):
        n_slices, n_particles = x.shape
        xij2 = zeros((n_slices,n_particles,n_particles),dtype=double,order='Fortran')
        vpi.xij.update_xij(xij2,x,1,n_slices)
        s = s % n_slices
        i = i % n_particles
        j = j % n_particles
        self.assert_(xij2[s,i,j] == norm(x[s,i]-x[s,j])**2)

if __name__ == "__main__":
    unittest.main()
#@-node:gcross.20090806151612.1843:@thin test.py
#@-leo
