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
#@+node:gcross.20090916114839.1816:compute_partial_sum
class compute_partial_sum(unittest.TestCase):
    #@    @+others
    #@+node:gcross.20090916114839.1817:test_correctnes
    @with_checker
    def test_correctness(self, n1 = irange(1,6), n2 = irange(1,8)):
        number_of_particles = max(n1,n2)
        number_excited = min(n1,n2)
        excluded_index = randint(0,number_of_particles-1)
        amplitudes = exp(rand(number_of_particles)+1j*rand(number_of_particles))
        computed_partial_sum = 0
        for selected_indices in combinations(range(number_of_particles),number_excited):
            if excluded_index in selected_indices:
                computed_partial_sum += prod(amplitudes[list(selected_indices)])
        formula_partial_sum = vpif.angular_momentum.compute_partial_sum(amplitudes,excluded_index+1,number_excited)
        self.assertAlmostEqual(computed_partial_sum,formula_partial_sum)
    #@-node:gcross.20090916114839.1817:test_correctnes
    #@-others
#@-node:gcross.20090916114839.1816:compute_partial_sum
#@+node:gcross.20090817102318.1730:accumulate_gradient_fancy
class accumulate_gradient_fancy(unittest.TestCase):
    #@    @+others
    #@+node:gcross.20090825141639.1536:Basic properties
    #@+node:gcross.20090817102318.1731:test_finite
    @with_checker
    def test_finite(self,
            n_slices = irange(1,4),
            n_particles = irange(1,10),
            fixed_rotation_axis = irange(1,3),
        ):
        x = rand(n_slices,n_particles,3)
        N_rotating_particles = randint(0,n_particles)
        rotation_plane_axis_1, rotation_plane_axis_2 = vpif.angular_momentum.get_rotation_plane_axes(fixed_rotation_axis)
        gradient_phase = zeros(x.shape,dtype=double,order='Fortran')
        vpif.angular_momentum.accumulate_gradient_fancy(
            x,
            N_rotating_particles,
            rotation_plane_axis_1,rotation_plane_axis_2,
            gradient_phase
            )

        self.assert_(isfinite(gradient_phase).all())
    #@-node:gcross.20090817102318.1731:test_finite
    #@-node:gcross.20090825141639.1536:Basic properties
    #@+node:gcross.20090825141639.1537:Correctness
    #@+node:gcross.20090819152718.1587:test_correct
    @with_checker(number_of_calls=10)
    def test_correct(self,
            N_particles  = irange(1,5),
            N_dimensions = irange(2,5),
        ):
        N_rotating_particles = randint(0,N_particles)
        x = rand(1,N_particles,N_dimensions)
        rotation_plane_axis_1 = randint(1,N_dimensions-1)
        rotation_plane_axis_2 = randint(rotation_plane_axis_1+1,N_dimensions)
        angles = arctan2(
            x[0,:,rotation_plane_axis_2-1],
            x[0,:,rotation_plane_axis_1-1]
        )
        gradient_phase = zeros(x.shape,dtype=double,order='Fortran')
        vpif.angular_momentum.accumulate_gradient_fancy(
            x,
            N_rotating_particles,
            rotation_plane_axis_1,rotation_plane_axis_2,
            gradient_phase
            )
        for i in xrange(len(angles)):
            numerical_derivative = derivative(
                self.make_phase1(N_rotating_particles,angles,i),
                angles[i],
                dx=1e-6,
                n=1,
                order=13
            )
            rot_x = x[0,i,rotation_plane_axis_1-1]
            rot_y = x[0,i,rotation_plane_axis_2-1]
            rot_r_squared = rot_x**2 + rot_y**2
            self.assertAlmostEqual(
               -numerical_derivative * rot_y/rot_r_squared,
                gradient_phase[0,i,rotation_plane_axis_1-1]
            )
            self.assertAlmostEqual(
               +numerical_derivative * rot_x/rot_r_squared,
                gradient_phase[0,i,rotation_plane_axis_2-1]
            )
    #@-node:gcross.20090819152718.1587:test_correct
    #@+node:gcross.20090819152718.1588:phase
    @staticmethod
    def phase(N_rotating_particles, angles):
        C = 0
        S = 0
        for a in combinations(angles,N_rotating_particles):
            C += cos(sum(a))
            S += sin(sum(a))
        if C == -0:
            return -pi/2
        if C == +0:
            return pi/2
        else:
            return atan(S/C)
    #@-node:gcross.20090819152718.1588:phase
    #@+node:gcross.20090819152718.1589:make_phase1
    def make_phase1(self,N_rotating_particles,angles,index):
        angles = angles.copy()
        def phase1(angle):
            angles[index] = angle
            return self.phase(N_rotating_particles,angles)
        return phase1
    #@-node:gcross.20090819152718.1589:make_phase1
    #@-node:gcross.20090825141639.1537:Correctness
    #@-others
#@-node:gcross.20090817102318.1730:accumulate_gradient_fancy
#@+node:gcross.20090915142144.1662:compute_gradient_fancy_amplitude
class compute_gradient_fancy_amplitude(unittest.TestCase):
    #@    @+others
    #@+node:gcross.20090915142144.1663:Basic properties
    #@+node:gcross.20090915142144.1664:test_finite
    @with_checker
    def test_finite(self,
            n_particles = irange(1,10),
        ):
        x = rand(n_particles,2)
        N_rotating_particles = randint(0,n_particles)
        self.assert_(isfinite(
            vpif.angular_momentum.compute_gradient_fancy_amplitude(
                x,
                N_rotating_particles,
                1,2
                )
        ).all())
    #@-node:gcross.20090915142144.1664:test_finite
    #@-node:gcross.20090915142144.1663:Basic properties
    #@+node:gcross.20090915142144.1665:Correctness
    #@+node:gcross.20090915142144.1666:test_correct
    @with_checker(number_of_calls=10)
    def test_correct(self,
            N_particles = irange(1,5),
        ):
        N_rotating_particles = randint(0,N_particles)
        x = rand(N_particles,2)
        angles = arctan2(
            x[:,1],
            x[:,0]
        )
        gradient_amplitude = \
            vpif.angular_momentum.compute_gradient_fancy_amplitude(
                x,
                N_rotating_particles,
                1,2
            )
        for i in xrange(N_particles):
            rho = norm(x[i])
            angles_copy = angles.copy()
            def C(angle):
                angles_copy[i] = angle
                return __builtin__.sum(imap(cos,imap(__builtin__.sum,combinations(angles_copy,N_rotating_particles))))
            numerical_derivative = derivative(
                C,
                angles[i],
                dx=1e-6,
                n=1,
                order=13
            )*rho
            self.assertAlmostEqual(
                numerical_derivative,
                gradient_amplitude[i,0]
            )
            def S(angle):
                angles_copy[i] = angle
                return __builtin__.sum(imap(sin,imap(__builtin__.sum,combinations(angles_copy,N_rotating_particles))))
            numerical_derivative = derivative(
                S,
                angles[i],
                dx=1e-6,
                n=1,
                order=13
            )*rho
            self.assertAlmostEqual(
                numerical_derivative,
                gradient_amplitude[i,1]
            )
    #@-node:gcross.20090915142144.1666:test_correct
    #@-node:gcross.20090915142144.1665:Correctness
    #@-others
#@-node:gcross.20090915142144.1662:compute_gradient_fancy_amplitude
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
    @with_checker(number_of_calls=10)
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
#@+node:gcross.20090813184545.1726:compute_rotational_potential
class accumulate_rotation_potential(unittest.TestCase):
    #@    @+others
    #@+node:gcross.20090817102318.1733:test_finite
    @with_checker
    def test_finite(self,
            lambda_ = unit_interval_float,
            n_particles = irange(1,10),
            fixed_rotation_axis = irange(1,3),
            frame_angular_velocity=unit_interval_float,
        ):
        n_dimensions = 3
        rotation_plane_axis_1, rotation_plane_axis_2 = vpif.angular_momentum.get_rotation_plane_axes(fixed_rotation_axis)
        x = rand(n_particles,n_dimensions)
        N_rotating_particles = randint(0,n_particles)
        U = zeros((n_particles,),dtype=double,order='Fortran')
        first_angular_derivatives = rand(n_particles)
        vpif.angular_momentum.accumulate_rotation_potential(
            x,lambda_,
            first_angular_derivatives,
            rotation_plane_axis_1,rotation_plane_axis_2,
            frame_angular_velocity,
            U
        )
        self.assert_(isfinite(U).all())
    #@nonl
    #@-node:gcross.20090817102318.1733:test_finite
    #@+node:gcross.20090817102318.1726:test_clumping
    @with_checker
    def test_clumping(self,
            n_particles=irange(2,10),
            angular_width_1=frange(0,0.6),
            angular_width_2=frange(0,0.6),
        ):
        n_dimensions = 3
        x = zeros((n_particles,n_dimensions),dtype=double,order='Fortran')
        U = zeros((n_particles),dtype=double,order='Fortran')
        N_rotating_particles = randint(1,n_particles)
        frame_angular_velocity = 0.0
        fixed_rotation_axis = 3
        rotation_plane_axis_1, rotation_plane_axis_2 = vpif.angular_momentum.get_rotation_plane_axes(fixed_rotation_axis)
        potentials = []
        for angular_width in sorted([angular_width_1,angular_width_2]):
            angles = (array(range(n_particles))*angular_width).reshape(1,n_particles)
            radii = rand(n_particles)
            x[:,0] = radii*cos(angles)
            x[:,1] = radii*sin(angles)
            x[:,2] = rand(n_particles)
            first_angular_derivatives = vpif.angular_momentum.compute_angular_derivatives(
                x,
                rotation_plane_axis_1, rotation_plane_axis_2, N_rotating_particles
            )
            U = zeros((n_particles,),dtype=double,order='Fortran')
            vpif.angular_momentum.accumulate_rotation_potential(
                x,0.5,
                first_angular_derivatives,
                rotation_plane_axis_1,rotation_plane_axis_2,
                frame_angular_velocity,
                U
            )
            potentials.append((angular_width,sum(U)))
        self.assert_(potentials[1] >= potentials[0])
    #@nonl
    #@-node:gcross.20090817102318.1726:test_clumping
    #@+node:gcross.20090817102318.1751:test_angular_momentum_raises_potential
    @with_checker
    def test_angular_momentum_raises_potential(self,
            n_particles = irange(1,10),
            fixed_rotation_axis = irange(1,3),
        ):
        n_dimensions = 3
        rotation_plane_axis_1, rotation_plane_axis_2 = vpif.angular_momentum.get_rotation_plane_axes(fixed_rotation_axis)
        x = rand(n_particles,n_dimensions)
        potentials = []
        for N_rotating_particles in [0,randint(1,n_particles)]:
            first_angular_derivatives = vpif.angular_momentum.compute_angular_derivatives(
                x,
                rotation_plane_axis_1, rotation_plane_axis_2, N_rotating_particles
            )
            U = zeros((n_particles,),dtype=double,order='Fortran')
            vpif.angular_momentum.accumulate_rotation_potential(
                x,0.5,
                first_angular_derivatives,
                rotation_plane_axis_1,rotation_plane_axis_2,
                0,
                U
            )
            potentials.append(sum(U))
        self.assert_(potentials[0] < potentials[1])
    #@nonl
    #@-node:gcross.20090817102318.1751:test_angular_momentum_raises_potential
    #@+node:gcross.20090817102318.1753:test_angular_momentum_cancels_frame_rotation
    @with_checker
    def test_angular_momentum_cancels_frame_rotation(self,
            n_particles = irange(1,10),
            fixed_rotation_axis = irange(1,3),
            frame_angular_velocity=unit_interval_float,
        ):
        n_dimensions = 3
        rotation_plane_axis_1, rotation_plane_axis_2 = vpif.angular_momentum.get_rotation_plane_axes(fixed_rotation_axis)
        x = rand(n_particles,3)
        for i in xrange(n_particles):
            x[i,fixed_rotation_axis-1] = 0
            x[i] /= norm(x[i])
        potentials = []
        preferred_momentum = round(frame_angular_velocity*n_particles)
        chosen_momentum = preferred_momentum
        while(chosen_momentum == preferred_momentum):
            chosen_momentum = randint(0,n_particles)
        for N_rotating_particles in [preferred_momentum,chosen_momentum]:
            first_angular_derivatives = vpif.angular_momentum.compute_angular_derivatives(
                x,
                rotation_plane_axis_1, rotation_plane_axis_2, N_rotating_particles
            )
            U = zeros((n_particles,),dtype=double,order='Fortran')
            vpif.angular_momentum.accumulate_rotation_potential(
                x,0.5,
                first_angular_derivatives,
                rotation_plane_axis_1,rotation_plane_axis_2,
                frame_angular_velocity,
                U
            )
            potentials.append(sum(U))
        self.assert_(potentials[0] < potentials[1])
    #@nonl
    #@-node:gcross.20090817102318.1753:test_angular_momentum_cancels_frame_rotation
    #@-others
#@-node:gcross.20090813184545.1726:compute_rotational_potential
#@-others

tests = [
    get_rotation_plane_axes,
    perform_special_matmul,
    sum_over_symmetrizations,
    accumulate_gradient_fancy,
    compute_gradient_fancy_amplitude,
    compute_amps_and_sum_syms,
    compute_partial_sum,
    estimate_distance_to_node,
    #accumulate_rotation_potential
    ]

if __name__ == "__main__":
    unittest.main()
#@nonl
#@-node:gcross.20090807144330.2253:@thin angular_momentum.py
#@-leo
