#@+leo-ver=4-thin
#@+node:gcross.20090807144330.1674:@thin __init__.py
import unittest
from paycheck import generator
from numpy import array, zeros, all, double
from numpy.linalg import norm
from numpy.random import rand
import random
import vpi

#@+others
#@+node:gcross.20090807144330.2153:ArrayGenerator
# ------------------------------------------------------------------------------
# NumPy Support
# ------------------------------------------------------------------------------

class ArrayGenerator(generator.PayCheckGenerator):
    def __init__(self, example):
        self.shape = example.shape
        self.dtype = example.dtype

    def next(self):
        return array(rand(*[random.randint(1, generator.LIST_LEN) for _ in self.shape]),dtype=self.dtype)

generator.container_generators[type(zeros(()))] = ArrayGenerator
#@nonl
#@-node:gcross.20090807144330.2153:ArrayGenerator
#@+node:gcross.20090807144330.2155:specialized type generators
particle_paths_type = zeros((0,0),dtype=double)
#@-node:gcross.20090807144330.2155:specialized type generators
#@+node:gcross.20090807144330.2156:runner
def run_tests():
    tests = []
    for module in modules:
        for test_case in module.tests:
            tests.append(unittest.defaultTestLoader.loadTestsFromTestCase(test_case))
    unittest.TextTestRunner(verbosity=2).run(unittest.TestSuite(tests))


if __name__ == "__main__":
    run_tests()
#@-node:gcross.20090807144330.2156:runner
#@+node:gcross.20090807144330.2154:tests
from tests import xij, thermalize, angular_momentum, gfn, sample

modules = [
    gfn,
    thermalize,
    xij,
    angular_momentum,
    sample,
    ]
#@nonl
#@-node:gcross.20090807144330.2154:tests
#@-others
#@-node:gcross.20090807144330.1674:@thin __init__.py
#@-leo
