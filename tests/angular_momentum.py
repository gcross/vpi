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
import vpif

#@+others
#@+node:gcross.20090807171924.1721:get_rotation_plane_axes
class get_rotation_plane_axes(unittest.TestCase):
    def test_X_axis(self):
        self.assertEqual(vpif.angular_momentum.get_rotation_plane_axes(1),(2,3))
    def test_Y_axis(self):
        self.assertEqual(vpif.angular_momentum.get_rotation_plane_axes(2),(1,3))
    def test_Z_axis(self):
        self.assertEqual(vpif.angular_momentum.get_rotation_plane_axes(3),(1,2))
#@nonl
#@-node:gcross.20090807171924.1721:get_rotation_plane_axes
#@+node:gcross.20090807144330.2254:perform_special_matmul
class perform_special_matmul(unittest.TestCase):
    #@    @+others
    #@+node:gcross.20090817102318.1728:test_trivial_case
    def test_trivial_case(self):
        vector = array([2+0j],dtype=complex128)
        amplitudes = array([1+1j],dtype=complex128)
        vpif.angular_momentum.perform_special_matmul(vector,amplitudes)
        self.assertEqual(vector,array([2+2j],dtype=complex128))
    #@nonl
    #@-node:gcross.20090817102318.1728:test_trivial_case
    #@+node:gcross.20090817102318.1729:test_random_cases
    @with_checker
    def test_random_cases(self,values=[(complex,complex)]):
        if(len(values) == 0): return
        vector, amplitudes = (array(values_,dtype=complex128) for values_ in zip(*values))
        original_vector = vector.copy()
        vpif.angular_momentum.perform_special_matmul(vector,amplitudes)
        partial_sum = 0
        for i in xrange(len(vector)):
            partial_sum += original_vector[i]
            self.assertEqual(vector[i],partial_sum*amplitudes[i])
    #@nonl
    #@-node:gcross.20090817102318.1729:test_random_cases
    #@-others

#@-node:gcross.20090807144330.2254:perform_special_matmul
#@+node:gcross.20090807171924.1722:sum_over_symmetrizations
class sum_over_symmetrizations(unittest.TestCase):
    #@    @+others
    #@+node:gcross.20090817102318.1727:test_correctnes
    @with_checker
    def test_correctness(self, n1 = irange(1,6), n2 = irange(1,8)):
        number_of_particles = max(n1,n2)
        number_excited = min(n1,n2)
        amplitudes = exp(rand(number_of_particles)+1j*rand(number_of_particles))
        result = vpif.angular_momentum.sum_over_symmetrizations(amplitudes,number_excited)
        correct = __builtin__.sum(imap(prod,combinations(amplitudes,number_excited)))
        self.assertAlmostEqual(result,correct)
    #@nonl
    #@-node:gcross.20090817102318.1727:test_correctnes
    #@-others
#@-node:gcross.20090807171924.1722:sum_over_symmetrizations
#@+node:gcross.20090916114839.1816:compute_sum_and_its_derivatives
class compute_partial_sum(unittest.TestCase):
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
#@+node:gcross.20090915142144.1679:compute_amps_and_sum_syms
class compute_amps_and_sum_syms(unittest.TestCase):
    #@    @+others
    #@+node:gcross.20090915142144.1680:Basic properties
    #@+node:gcross.20090915142144.1681:test_finite
    @with_checker
    def test_finite(self,
            n_particles = irange(1,10),
        ):
        x = array(rand(n_particles,2),dtype=double,order='Fortran')
        n_rotating_particles = randint(0,n_particles)
        self.assert_(isfinite(
            vpif.angular_momentum.compute_amps_and_sum_syms(
                x,
                n_rotating_particles,
                1,2
                )
        ).all())
    #@-node:gcross.20090915142144.1681:test_finite
    #@-node:gcross.20090915142144.1680:Basic properties
    #@+node:gcross.20090915142144.1682:Correctness
    #@+node:gcross.20090915142144.1683:test_correct
    @with_checker
    def test_correct(self,
            n_particles = irange(1,10),
        ):
        x = array(rand(n_particles,2),dtype=double,order='Fortran')
        n_rotating_particles = randint(0,n_particles)
        computed_value = \
            vpif.angular_momentum.compute_amps_and_sum_syms(
                x,
                n_rotating_particles,
                1,2
                )
        C = 0
        S = 0
        for a in combinations(arctan2(x[:,1],x[:,0]),n_rotating_particles):
            C += cos(sum(a))
            S += sin(sum(a))
        correct_value = C+1j*S

        self.assertAlmostEqual(correct_value,computed_value)
    #@-node:gcross.20090915142144.1683:test_correct
    #@-node:gcross.20090915142144.1682:Correctness
    #@-others
#@-node:gcross.20090915142144.1679:compute_amps_and_sum_syms
#@+node:gcross.20090916153857.1676:estimate_distance_to_node
class estimate_distance_to_node(unittest.TestCase):
    #@    @+others
    #@+node:gcross.20090916153857.1677:Basic properties
    #@+node:gcross.20090916153857.1678:test_finite
    @with_checker
    def dont_test_finite(self,
            n_particles = irange(2,10),
        ):
        x = array(rand(n_particles,2),dtype=double,order='Fortran')
        n_rotating_particles = randint(0,n_particles)
        self.assert_(isfinite(
            vpif.angular_momentum.estimate_distance_to_node(
                x,
                n_rotating_particles,
                1,2
            )
        ).all())
    #@-node:gcross.20090916153857.1678:test_finite
    #@-node:gcross.20090916153857.1677:Basic properties
    #@+node:gcross.20090916153857.1679:Correctness
    #@+node:gcross.20090916153857.1680:test_correct
    @with_checker
    def test_correct(self,
            n_particles = irange(2,10),
        ):
        x = array(rand(n_particles,2),dtype=double,order='Fortran')
        n_rotating_particles = randint(0,n_particles)

        amplitude = vpif.angular_momentum.compute_amps_and_sum_syms(
            x,
            n_rotating_particles,
            1,2
        )

        correct_value = norm(lstsq(
            vpif.angular_momentum.compute_gradient_fancy_amplitude(
                x,
                n_rotating_particles,
                1,2
            ).transpose(),
            array([real(amplitude),imag(amplitude)]),
            cond=-1,
            overwrite_a=True,
            overwrite_b=True
        )[0])

        computed_value = \
            vpif.angular_momentum.estimate_distance_to_node(
                x,
                n_rotating_particles,
                1,2
            )

        self.assertAlmostEqual(correct_value,computed_value)
    #@-node:gcross.20090916153857.1680:test_correct
    #@-node:gcross.20090916153857.1679:Correctness
    #@-others
#@-node:gcross.20090916153857.1676:estimate_distance_to_node
#@-others

tests = [
    get_rotation_plane_axes,
    perform_special_matmul,
    sum_over_symmetrizations,
    compute_amps_and_sum_syms,
    compute_partial_sum,
    estimate_distance_to_node,
    ]

if __name__ == "__main__":
    unittest.main()
#@nonl
#@-node:gcross.20090807144330.2253:@thin angular_momentum.py
#@-leo
