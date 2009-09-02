#@+leo-ver=4-thin
#@+node:gcross.20090812093015.1742:@thin gfn.py
#@@language python
#@@tabwidth -4

import unittest
from paycheck import *
from numpy import zeros, double, dot
from numpy.linalg import norm
from numpy.random import rand
from random import randint
import vpif

#@+others
#@+node:gcross.20090812093015.1743:gfn2_sp
class gfn2_sp(unittest.TestCase):

    @with_checker
    def test_correctness(self, n_slices = irange(1,10), n_particles = irange(1,5), dt = float):
        sl_start = randint(1,n_slices)
        sl_end = randint(sl_start,n_slices)
        ip = randint(1,n_particles)
        U = rand(n_slices,n_particles)
        ln_gfn = vpif.gfn.gfn2_sp(sl_start,sl_end,ip,U,dt)
        self.assertAlmostEqual(ln_gfn,-dt*sum(U[sl_start-1:sl_end,ip-1]))
#@nonl
#@-node:gcross.20090812093015.1743:gfn2_sp
#@+node:gcross.20090812093015.1747:gfn4_sp
class gfn4_sp(unittest.TestCase):

    #@    @+others
    #@+node:gcross.20090812093015.1749:initialize_weights
    @staticmethod
    def initialize_weights(n_slices):
        assert (n_slices%2==0)
        c_slice = n_slices/2
        assert (c_slice%2==1)
        U_weight = zeros((n_slices+2))
        gU2_weight = zeros((n_slices+2))
        for k in xrange(1,c_slice+1):
            U_weight[k] = (k+1)%2+1
            gU2_weight[k] = (k+1)%2
        for k in xrange(c_slice+2,n_slices+1):
            U_weight[k] = k%2+1
            gU2_weight[k] = k%2
        U_weight[1] = 0.5
        U_weight[n_slices] = 0.5
        U_weight[c_slice] = 0.5
        U_weight[c_slice+1] = 0.5
        U_weight = U_weight[1:n_slices+1]
        gU2_weight = gU2_weight[1:n_slices+1]

        return (U_weight, gU2_weight)
    #@-node:gcross.20090812093015.1749:initialize_weights
    #@+node:gcross.20090812093015.1752:test_correctness
    @with_checker
    def test_correctness(self,c_slice = irange(1,15,2),n_particles = irange(1,5),lam = positive_float(1e-5,100),dt = positive_float(1e-5,100)):
        n_slices = c_slice*2
        self.assertEqual(1,c_slice%2)

        (U_weight, gU2_weight) = self.initialize_weights(n_slices)

        sl_start = randint(1,n_slices)
        sl_end = randint(sl_start,n_slices)
        ip = randint(1,n_particles)
        U = rand(n_slices,n_particles)
        gradU2 = rand(n_slices)
        ln_gfn = vpif.gfn.gfn4_sp(sl_start,sl_end,ip,U,gradU2,U_weight,gU2_weight,lam,dt)

        self.assertAlmostEqual(ln_gfn,
            -2.0*dt*dot(
                U[sl_start-1:sl_end,ip-1],
                U_weight[sl_start-1:sl_end]
            )/3.0
            -2.0*lam*(dt**3)*dot(
                gradU2[sl_start-1:sl_end],
                gU2_weight[sl_start-1:sl_end]
            )/9.0
        )
    #@nonl
    #@-node:gcross.20090812093015.1752:test_correctness
    #@-others
#@-node:gcross.20090812093015.1747:gfn4_sp
#@+node:gcross.20090828095451.1459:gfn_hard_wall_contribution
class gfn_hard_wall_contribution(unittest.TestCase):

    @with_checker
    def test_finite(self,
            n_slices = irange(1,10),
            n_particles = irange(1,5),
            n_dimensions = irange(1,5),
            lambda_ = float,
            dt = float
        ):
        q = rand(n_slices,n_particles,n_dimensions)
        hard_wall_locations = rand(n_dimensions)
        slice_start = randint(1,n_slices)
        slice_end = randint(slice_start,n_slices)
        particle_number = randint(1,n_particles)
        vpif.gfn.gfn_hard_wall_contribution(
            q,
            hard_wall_locations,
            lambda_, dt,
            slice_start, slice_end,
            particle_number
        )
#@nonl
#@-node:gcross.20090828095451.1459:gfn_hard_wall_contribution
#@+node:gcross.20090828095451.1461:gfn_hard_sphere_contribution
class gfn_hard_sphere_contribution(unittest.TestCase):

    @with_checker
    def test_finite(self,
            n_slices = irange(1,10),
            n_particles = irange(1,5),
            dt = float,
            hard_sphere_radius = float,
        ):
        xij2 = rand(n_slices,n_particles,n_particles)
        slice_start = randint(1,n_slices)
        slice_end = randint(slice_start,n_slices)
        particle_number = randint(1,n_particles)
        vpif.gfn.gfn_hard_sphere_contribution(
            xij2,
            dt,hard_sphere_radius,
            slice_start, slice_end,
            particle_number
        )
#@nonl
#@-node:gcross.20090828095451.1461:gfn_hard_sphere_contribution
#@+node:gcross.20090828095451.1463:gfn_hard_sphere_contribution2
class gfn_hard_sphere_contribution2(unittest.TestCase):

    @with_checker
    def test_finite(self,
            n_slices = irange(1,10),
            n_particles = irange(1,5),
            n_dimensions = irange(1,5),
            lambda_ = float,
            dt = float,
            hard_sphere_radius = float,
        ):
        q = rand(n_slices,n_particles,n_dimensions)
        slice_start = randint(1,n_slices)
        slice_end = randint(slice_start,n_slices)
        particle_number = randint(1,n_particles)
        vpif.gfn.gfn_hard_sphere_contribution2(
            q,
            lambda_,dt,hard_sphere_radius,
            slice_start, slice_end,
            particle_number
        )
#@nonl
#@-node:gcross.20090828095451.1463:gfn_hard_sphere_contribution2
#@-others

tests = [
    gfn2_sp,
    gfn4_sp,
    gfn_hard_wall_contribution,
    gfn_hard_sphere_contribution,
    gfn_hard_sphere_contribution2,
    ]

if __name__ == "__main__":
    unittest.main()
#@nonl
#@-node:gcross.20090812093015.1742:@thin gfn.py
#@-leo
