#@+leo-ver=4-thin
#@+node:gcross.20090807144330.1674:@thin __init__.py
import unittest
from paycheck import generator
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
        return array(rand(*[random.randint(1, generator.LIST_LEN) for _ in self.shape]),dtype=self.dtype)

generator.container_generators[type(zeros(()))] = ArrayGenerator

particle_paths_type = zeros((0,0),dtype=double)

from tests import xij

modules = [
    xij,
    ]

def run_tests():
    tests = []
    for module in modules:
        for test_case in module.tests:
            tests.append(unittest.defaultTestLoader.loadTestsFromTestCase(test_case))
    unittest.TextTestRunner(verbosity=2).run(unittest.TestSuite(tests))


if __name__ == "__main__":
    run_tests()
#@-node:gcross.20090807144330.1674:@thin __init__.py
#@-leo
