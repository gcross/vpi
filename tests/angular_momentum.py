#@+leo-ver=4-thin
#@+node:gcross.20090807144330.2253:@thin angular_momentum.py
#@@language python
#@@tabwidth -4

import unittest
from paycheck import *
from numpy import *
from numpy.random import rand
from numpy.linalg import norm
from random import randint
from itertools import imap, combinations
import __builtin__
import vpi

#@+others
#@+node:gcross.20090807171924.1721:get_rotation_plane_axes
class get_rotation_plane_axes(unittest.TestCase):
    def test_X_axis(self):
        self.assertEqual(vpi.angular_momentum.get_rotation_plane_axes(1),(2,3))
    def test_Y_axis(self):
        self.assertEqual(vpi.angular_momentum.get_rotation_plane_axes(2),(1,3))
    def test_Z_axis(self):
        self.assertEqual(vpi.angular_momentum.get_rotation_plane_axes(3),(1,2))
#@-node:gcross.20090807171924.1721:get_rotation_plane_axes
#@+node:gcross.20090807144330.2254:perform_special_matmul
class perform_special_matmul(unittest.TestCase):
    #@    @+others
    #@+node:gcross.20090817102318.1728:test_trivial_case
    def test_trivial_case(self):
        vector = array([2+0j],dtype=complex128)
        amplitudes = array([1+1j],dtype=complex128)
        vpi.angular_momentum.perform_special_matmul(vector,amplitudes)
        self.assertEqual(vector,array([2+2j],dtype=complex128))
    #@-node:gcross.20090817102318.1728:test_trivial_case
    #@+node:gcross.20090817102318.1729:test_random_cases
    @with_checker
    def test_random_cases(self,values=[(complex,complex)]):
        if(len(values) == 0): return
        vector, amplitudes = (array(values_,dtype=complex128) for values_ in zip(*values))
        original_vector = vector.copy()
        vpi.angular_momentum.perform_special_matmul(vector,amplitudes)
        partial_sum = 0
        for i in xrange(len(vector)):
            partial_sum += original_vector[i]
            self.assertEqual(vector[i],partial_sum*amplitudes[i])
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
        result = vpi.angular_momentum.sum_over_symmetrizations(amplitudes,number_excited)
        correct = __builtin__.sum(imap(prod,combinations(amplitudes,number_excited)))
        self.assertAlmostEqual(result,correct)
    #@-node:gcross.20090817102318.1727:test_correctnes
    #@-others
#@-node:gcross.20090807171924.1722:sum_over_symmetrizations
#@+node:gcross.20090817102318.1730:compute_angular_derivatives
class compute_angular_derivatives(unittest.TestCase):
    #@    @+others
    #@+node:gcross.20090817102318.1731:test_finite
    @with_checker
    def test_finite(self,
            n_particles = irange(1,10),
            fixed_rotation_axis = irange(1,3),
        ):
        x = rand(n_particles,3)
        N_rotating_particles = randint(0,n_particles)
        rotation_plane_axis_1, rotation_plane_axis_2 = vpi.angular_momentum.get_rotation_plane_axes(fixed_rotation_axis)
        first_derivatives, second_derivatives = vpi.angular_momentum.compute_angular_derivatives(x,rotation_plane_axis_1,rotation_plane_axis_2,N_rotating_particles)
        self.assert_(isfinite(first_derivatives).all())
        self.assert_(isfinite(second_derivatives).all())
    #@-node:gcross.20090817102318.1731:test_finite
    #@+node:gcross.20090818081913.1338:test_2nd_derivatives_symmetric
    @with_checker
    def test_2nd_derivatives_symmetric(self,
            n_particles = irange(1,10),
            fixed_rotation_axis = irange(1,3),
        ):
        x = rand(n_particles,3)
        N_rotating_particles = randint(0,n_particles)
        rotation_plane_axis_1, rotation_plane_axis_2 = vpi.angular_momentum.get_rotation_plane_axes(fixed_rotation_axis)
        first_derivatives, second_derivatives = vpi.angular_momentum.compute_angular_derivatives(x,rotation_plane_axis_1,rotation_plane_axis_2,N_rotating_particles)
        for i in xrange(n_particles):
            for j in xrange(n_particles):
                self.assertAlmostEqual(second_derivatives[i,j],second_derivatives[j,i])
    #@-node:gcross.20090818081913.1338:test_2nd_derivatives_symmetric
    #@-others
#@-node:gcross.20090817102318.1730:compute_angular_derivatives
#@+node:gcross.20090813184545.1726:compute_effective_rotational_potential
class compute_effective_rotational_potential(unittest.TestCase):
    #@    @+others
    #@+node:gcross.20090817102318.1733:test_finite
    @with_checker
    def test_finite(self,
            n_slices = irange(1,5),
            n_particles = irange(1,10),
            fixed_rotation_axis = irange(1,3),
            frame_angular_velocity=unit_interval_float,
        ):
        n_dimensions = 3
        rotation_plane_axis_1, rotation_plane_axis_2 = vpi.angular_momentum.get_rotation_plane_axes(fixed_rotation_axis)
        x = rand(n_slices,n_particles,n_dimensions)
        N_rotating_particles = randint(0,n_particles)
        U = zeros((n_slices,n_particles),dtype=double,order='Fortran')
        gradU = zeros((n_slices,n_particles,n_dimensions),dtype=double,order='Fortran')
        vpi.angular_momentum. compute_effective_rotational_potential(x,rotation_plane_axis_1,rotation_plane_axis_2,frame_angular_velocity,N_rotating_particles,U,gradU)
        self.assert_(isfinite(U).all())
    #@nonl
    #@-node:gcross.20090817102318.1733:test_finite
    #@+node:gcross.20090817102318.1726:test_clumping
    @with_checker
    def test_clumping(self,
            n_slices=irange(1,4),
            n_particles=irange(2,10),
            angular_width_1=frange(0,0.6),
            angular_width_2=frange(0,0.6),
        ):
        n_dimensions = 3
        x = zeros((n_slices,n_particles,n_dimensions),dtype=double,order='Fortran')
        U = zeros((n_slices,n_particles),dtype=double,order='Fortran')
        gradU = zeros((n_slices,n_particles,n_dimensions),dtype=double,order='Fortran')
        N_rotating_particles = randint(1,n_particles)
        frame_angular_velocity = 0.0
        fixed_rotation_axis = 3
        rotation_plane_axis_1, rotation_plane_axis_2 = vpi.angular_momentum.get_rotation_plane_axes(fixed_rotation_axis)
        potentials = []
        for angular_width in sorted([angular_width_1,angular_width_2]):
            angles = (array(range(n_particles))*angular_width).reshape(1,n_particles)
            radii = rand(n_slices,n_particles)
            x[:,:,0] = radii*cos(angles)
            x[:,:,1] = radii*sin(angles)
            x[:,:,2] = rand(n_slices,n_particles)
            U[...] = 0
            gradU[...] = 0
            vpi.angular_momentum.compute_effective_rotational_potential(
                x,
                rotation_plane_axis_1,rotation_plane_axis_2,frame_angular_velocity,N_rotating_particles,
                U,gradU
            )
            potentials.append((angular_width,sum(U)))
        self.assert_(potentials[1] >= potentials[0])
    #@nonl
    #@-node:gcross.20090817102318.1726:test_clumping
    #@+node:gcross.20090817102318.1751:test_angular_momentum_raises_potential
    @with_checker
    def test_angular_momentum_raises_potential(self,
            n_slices = irange(1,5),
            n_particles = irange(1,10),
            fixed_rotation_axis = irange(1,3),
        ):
        n_dimensions = 3
        rotation_plane_axis_1, rotation_plane_axis_2 = vpi.angular_momentum.get_rotation_plane_axes(fixed_rotation_axis)
        x = rand(n_slices,n_particles,n_dimensions)
        potentials = []
        for N_rotating_particles in [0,randint(1,n_particles)]:
            U = zeros((n_slices,n_particles),dtype=double,order='Fortran')
            gradU = zeros((n_slices,n_particles,n_dimensions),dtype=double,order='Fortran')
            vpi.angular_momentum.compute_effective_rotational_potential(x,rotation_plane_axis_1,rotation_plane_axis_2,0,N_rotating_particles,U,gradU)
            potentials.append(sum(U))
        self.assert_(potentials[0] < potentials[1])
    #@nonl
    #@-node:gcross.20090817102318.1751:test_angular_momentum_raises_potential
    #@+node:gcross.20090817102318.1753:test_angular_momentum_cancels_frame_rotation
    @with_checker
    def test_angular_momentum_cancels_frame_rotation(self,
            n_slices = irange(1,5),
            n_particles = irange(1,10),
            fixed_rotation_axis = irange(1,3),
            frame_angular_velocity=unit_interval_float,
        ):
        n_dimensions = 3
        rotation_plane_axis_1, rotation_plane_axis_2 = vpi.angular_momentum.get_rotation_plane_axes(fixed_rotation_axis)
        x = rand(n_slices,n_particles,3)
        for i in xrange(n_slices):
            for j in xrange(n_particles):
                x[i,j,fixed_rotation_axis-1] = 0
                x[i,j] /= norm(x[i,j])
        potentials = []
        preferred_momentum = round(frame_angular_velocity*n_particles)
        chosen_momentum = preferred_momentum
        while(chosen_momentum == preferred_momentum):
            chosen_momentum = randint(0,n_particles)
        for N_rotating_particles in [preferred_momentum,chosen_momentum]:
            U = zeros((n_slices,n_particles),dtype=double,order='Fortran')
            gradU = zeros((n_slices,n_particles,n_dimensions),dtype=double,order='Fortran')
            vpi.angular_momentum.compute_effective_rotational_potential(x,rotation_plane_axis_1,rotation_plane_axis_2,frame_angular_velocity,N_rotating_particles,U,gradU)
            potentials.append(sum(U))
        self.assert_(potentials[0] < potentials[1])
    #@nonl
    #@-node:gcross.20090817102318.1753:test_angular_momentum_cancels_frame_rotation
    #@-others
#@-node:gcross.20090813184545.1726:compute_effective_rotational_potential
#@-others

tests = [
    get_rotation_plane_axes,
    perform_special_matmul,
    sum_over_symmetrizations,
    compute_angular_derivatives,
    compute_effective_rotational_potential
    ]

if __name__ == "__main__":
    unittest.main()
#@-node:gcross.20090807144330.2253:@thin angular_momentum.py
#@-leo
