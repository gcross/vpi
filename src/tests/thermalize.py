#@+leo-ver=4-thin
#@+node:gcross.20090807144330.2151:@thin thermalize.py
#@@language python
#@@tabwidth -4

import unittest
from paycheck import with_checker
from paycheck.generator import positive_float, non_negative_float, irange
from numpy import array, zeros, double, float64
from numpy.linalg import norm
from numpy.random import rand
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
            irange(1,5),
        number_of_calls=10
    )
    def test_null_case(self,
            n_slices,
            n_particles,
            n_dimensions,
            n_dimensions_rotation
        ):
        x = array(rand(n_slices,n_particles,n_dimensions),order='Fortran')
        xij2 = array(rand(n_slices,n_particles,n_particles),order='Fortran')
        x_rot = array(zeros((n_slices,n_particles,n_dimensions_rotation)),order='Fortran')
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
            x_rot=x_rot,
            Uij_rot_func=null_func
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
            irange(1,5),
        number_of_calls=10
    )
    def test_constant_case(self,
            n_slices,
            n_particles,
            n_dimensions,
            n_dimensions_rotation
        ):
        x = array(rand(n_slices,n_particles,n_dimensions),order='Fortran')
        xij2 = array(rand(n_slices,n_particles,n_particles),order='Fortran')
        x_rot = array(zeros((n_slices,n_particles,n_dimensions_rotation)),order='Fortran')
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
            x_rot=x_rot,
            Uij_rot_func=null_func
            )
        self.failIf(reject_flag)
        self.assert_((U==2).all())
        self.assert_((gradU2==0).all())
    #@-node:gcross.20090812093015.1721:test_constant_case
    #@-others
#@-node:gcross.20090812093015.1718:compute_potential
#@-others

tests = [
    accept_path,
    compute_physical_potential,
    ]

if __name__ == "__main__":
    unittest.main()
#@-node:gcross.20090807144330.2151:@thin thermalize.py
#@-leo
