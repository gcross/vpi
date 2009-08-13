#@+leo-ver=4-thin
#@+node:gcross.20090807144330.2151:@thin thermalize.py
#@@language python
#@@tabwidth -4

import unittest
from paycheck import with_checker
from paycheck.generator import positive_float, non_negative_float, irange, frange
from numpy import array, zeros, double, float64, isfinite
from numpy.linalg import norm
from numpy.random import rand
from random import randint
from tests import particle_paths_type
import vpi

#@+others
#@+node:gcross.20090807144330.2152:accept_path
class accept_path(unittest.TestCase):
    @with_checker(float,float)
    def test_always_accepts_minimizing_move(self,p1,p2):
        self.assert_(
            (vpi.thermalize.accept_path(p1,p2) == 0) or
            (p1 < p2)
        )
#@-node:gcross.20090807144330.2152:accept_path
#@+node:gcross.20090812093015.1718:compute_potential
class compute_physical_potential(unittest.TestCase):
    #@    @+others
    #@+node:gcross.20090812093015.1719:test_null_case
    @with_checker(
            irange(1,5),
            irange(1,5),
            irange(1,5),
        number_of_calls=10
    )
    def test_null_case(self,
            n_slices,
            n_particles,
            n_dimensions,
        ):
        x = array(rand(n_slices,n_particles,n_dimensions),order='Fortran')
        xij2 = array(rand(n_slices,n_particles,n_particles),order='Fortran')
        def null_func(*args):
            return double(0.0)
        def Uij_func(*args):
            return 0.0, False
        gnull = array(zeros((n_particles,n_dimensions)),dtype=double,order='Fortran')
        def null_grad_func(*args):
            return gnull
        U,gradU2,reject_flag = vpi.thermalize.compute_physical_potential(
            x,xij2,
            null_func,
            null_grad_func,
            Uij_func,
            null_grad_func,
            1,n_slices,
            )
        self.failIf(reject_flag)
        self.assert_((U==0).all())
        self.assert_((gradU2==0).all())
    #@nonl
    #@-node:gcross.20090812093015.1719:test_null_case
    #@+node:gcross.20090812093015.1721:test_constant_case
    @with_checker(
            irange(1,5),
            irange(1,5),
            irange(1,5),
        number_of_calls=10
    )
    def test_constant_case(self,
            n_slices,
            n_particles,
            n_dimensions,
        ):
        x = array(rand(n_slices,n_particles,n_dimensions),order='Fortran')
        xij2 = array(rand(n_slices,n_particles,n_particles),order='Fortran')
        def null_func(*args):
            return double(0.0)
        def constant_func(*args):
            return double(1.0)
        def Uij_func(*args):
            return double(1.0), False
        gnull = array(zeros((n_particles,n_dimensions)),dtype=double,order='Fortran')
        def null_grad_func(*args):
            return gnull
        U,gradU2,reject_flag = vpi.thermalize.compute_physical_potential(
            x,xij2,
            constant_func,
            null_grad_func,
            Uij_func,
            null_grad_func,
            1,n_slices,
            )
        self.failIf(reject_flag)
        self.assert_((U==2).all())
        self.assert_((gradU2==0).all())
    #@-node:gcross.20090812093015.1721:test_constant_case
    #@-others
#@-node:gcross.20090812093015.1718:compute_potential
#@+node:gcross.20090813093326.1730:compute_log_acceptance_weight
class compute_log_acceptance_weight(unittest.TestCase):
    #@    @+others
    #@+node:gcross.20090813095726.1739:test_null_case
    @with_checker(
            irange(1,11,2),
            irange(1,5),
            irange(1,5),
            irange(1,3),
            frange(0,1),
            frange(0,1),
            bool,
        )
    def test_null_case(self,
            c_slice,
            n_particles,
            n_dimensions,
            fixed_rotation_axis,
            lam,
            dtau,
            use_4th_order_green_function,
        ):
        n_slices = c_slice * 2
        x = array(rand(n_slices,n_particles,n_dimensions),order='Fortran')
        xij2 = array(rand(n_slices,n_particles,n_particles),order='Fortran')
        move_start = randint(1,n_slices)
        move_end = randint(move_start,n_slices)
        particle_number = randint(1,n_particles)
        def null_func(*args):
            return double(0.0)
        def Uij_func(*args):
            return 0.0, False
        gnull = array(zeros((n_particles,n_dimensions)),dtype=double,order='Fortran')
        def null_grad_func(*args):
            return gnull
        U_weights, gU2_weights = vpi.gfn.initialize_4th_order_weights(n_slices)
        U, gradU2, reject_flag, weight = vpi.thermalize.compute_log_acceptance_weight(
            x,xij2,
            move_start,move_end,
            particle_number,
            null_func,
            null_grad_func,
            Uij_func,
            null_grad_func,
            U_weights,gU2_weights,
            fixed_rotation_axis,
            0,0,
            lam, dtau,
            use_4th_order_green_function,
            null_func,
            null_func,
        )
        self.failIf(reject_flag)
        self.assert_((U==0).all())
        self.assert_((gradU2==0).all())
        self.assertEqual(0,weight)
    #@-node:gcross.20090813095726.1739:test_null_case
    #@+node:gcross.20090813095726.1741:test_that_angular_momentum_lowers_weight
    @with_checker(
            irange(3,11,2),
            irange(1,10),
            irange(1,3),
            frange(0,1),
            frange(0,1),
            bool,
        )
    def test_that_angular_momentum_lowers_weight(self,
            c_slice,
            n_particles,
            fixed_rotation_axis,
            lam,
            dtau,
            use_4th_order_green_function,
        ):
        fixed_angular_momentum = randint(1,n_particles)
        n_slices = c_slice * 2
        n_dimensions = 3
        x = array(rand(n_slices,n_particles,n_dimensions),order='Fortran')
        xij2 = array(rand(n_slices,n_particles,n_particles),order='Fortran')
        move_start = randint(1,n_slices)
        move_end = randint(move_start,n_slices)
        particle_number = randint(1,n_particles)
        def null_func(*args):
            return double(0.0)
        def Uij_func(*args):
            return 0.0, False
        gnull = array(zeros((n_particles,n_dimensions)),dtype=double,order='Fortran')
        def null_grad_func(*args):
            return gnull
        U_weights, gU2_weights = vpi.gfn.initialize_4th_order_weights(n_slices)
        U, gradU2, reject_flag, weight = vpi.thermalize.compute_log_acceptance_weight(
            x,xij2,
            move_start,move_end,
            particle_number,
            null_func,
            null_grad_func,
            Uij_func,
            null_grad_func,
            U_weights,gU2_weights,
            fixed_rotation_axis,0,fixed_angular_momentum,
            lam, dtau,
            use_4th_order_green_function,
            null_func,
            null_func,
        )
        self.failIf(reject_flag)
        U = U[move_start-1:move_end,particle_number-1]
        self.assert_(isfinite(U).all())
        self.assert_((U>0).all())
        self.assert_(isfinite(weight))
        self.assert_(weight<0)
    #@-node:gcross.20090813095726.1741:test_that_angular_momentum_lowers_weight
    #@-others
#@-node:gcross.20090813093326.1730:compute_log_acceptance_weight
#@-others

tests = [
    accept_path,
    compute_physical_potential,
    compute_log_acceptance_weight
    ]

if __name__ == "__main__":
    unittest.main()
#@-node:gcross.20090807144330.2151:@thin thermalize.py
#@-leo
