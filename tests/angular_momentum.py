#@+leo-ver=4-thin
#@+node:gcross.20090807144330.2253:@thin angular_momentum.py
#@@language python
#@@tabwidth -4

import unittest
from paycheck import *
from numpy import *
from numpy.random import rand
from numpy.linalg import norm
from scipy.linalg import lstsq
from scipy.misc import derivative
from random import randint
import itertools
from itertools import imap, combinations
from functools import partial
from math import atan
import __builtin__
import vpi.fortran as vpif

#@+others
#@+node:gcross.20090916114839.1816:compute_sum_and_its_derivatives
class compute_sum_and_its_derivatives(unittest.TestCase):
    #@    @+others
    #@+node:gcross.20090916114839.1817:test_correctnes
    @with_checker
    def test_correctness(self, number_of_amplitudes_to_include = irange(0,6), n = irange(1,6)):
        number_of_amplitudes = number_of_amplitudes_to_include+n
        amplitudes = rand(number_of_amplitudes) + 1j*rand(number_of_amplitudes)
        sum_over_symmetrizations, gradient_of_sum = vpif.angular_momentum.compute_sum_and_its_derivatives(amplitudes,number_of_amplitudes_to_include)
        correct_sum_over_symmetrizations = __builtin__.sum(product(seq) for seq in itertools.combinations(list(amplitudes),number_of_amplitudes_to_include))
        self.assertAlmostEqual(correct_sum_over_symmetrizations,sum_over_symmetrizations)
    #@-node:gcross.20090916114839.1817:test_correctnes
    #@-others
#@-node:gcross.20090916114839.1816:compute_sum_and_its_derivatives
#@-others

tests = [
    compute_sum_and_its_derivatives,
    ]

if __name__ == "__main__":
    unittest.main()
#@nonl
#@-node:gcross.20090807144330.2253:@thin angular_momentum.py
#@-leo
