#@+leo-ver=4-thin
#@+node:gcross.20090812093015.1742:@thin gfn.py
#@@language python
#@@tabwidth -4

import unittest
from paycheck import *
from numpy import *
from numpy.linalg import norm
from numpy.random import rand
from random import randint
import vpif

#@+others
#@+node:gcross.20090812093015.1743:gfn2_sp
class gfn2_sp(unittest.TestCase):
    #@    @+others
    #@+node:gcross.20100117204224.1702:test_correctness
    @with_checker
    def test_correctness(self, n_slices = irange(1,10), n_particles = irange(1,5), dt = float):
        sl_start = randint(1,n_slices)
        sl_end = randint(sl_start,n_slices)
        ip = randint(1,n_particles)
        U = rand(n_slices,n_particles)
        U_weights = rand(n_slices)
        ln_gfn = vpif.gfn.gfn2_sp(sl_start,sl_end,U,U_weights,dt)
        U_summed = sum(U,axis=1)
        self.assertAlmostEqual(ln_gfn,-dt*dot(U_summed[sl_start-1:sl_end],U_weights[sl_start-1:sl_end]))
    #@-node:gcross.20100117204224.1702:test_correctness
    #@-others
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
        ln_gfn = vpif.gfn.gfn4_sp(sl_start,sl_end,U,gradU2,U_weight,gU2_weight,lam,dt)
        U = sum(U,axis=1)

        self.assertAlmostEqual(ln_gfn,
            -2.0*dt*dot(
                U[sl_start-1:sl_end],
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
#@+node:gcross.20090918120753.2604:compute_green_fn_from_distances
#@@language python

class compute_green_fn_from_distances(unittest.TestCase):
    #@    @+others
    #@+node:gcross.20090918120753.2610:test_finite
    @with_checker
    def test_finite(self,
            n_slices = irange(6,100),
            magnitude = float,
            denominator = unit_interval_float,
        ):
        start_slice = randint(1,n_slices-1)
        end_slice = randint(start_slice+1,n_slices)
        distances = rand(n_slices) * magnitude
        self.assert_(isfinite(
            vpif.gfn.compute_green_fn_from_distances(
                distances,
                denominator,
                start_slice, end_slice,
            )
        ))
    #@-node:gcross.20090918120753.2610:test_finite
    #@+node:gcross.20090918120753.2616:test_finite_log
    @with_checker
    def test_finite_log(self,
            n_slices = irange(6,100),
            denominator = unit_interval_float,
        ):
        start_slice = randint(1,n_slices-1)
        end_slice = randint(start_slice+1,n_slices)
        distances = rand(n_slices)
        self.assert_(isfinite(log(
            vpif.gfn.compute_green_fn_from_distances(
                distances,
                denominator,
                start_slice, end_slice,
            )
        )))
    #@-node:gcross.20090918120753.2616:test_finite_log
    #@+node:gcross.20090918120753.2618:test_finite_log_half_path
    @with_checker
    def test_finite_log_half_path(self,
            n_slices = irange(6,100),
            denominator = unit_interval_float,
        ):
        start_slice = randint(n_slices/2+2,n_slices-1)
        end_slice = randint(start_slice+1,n_slices)
        distances = rand(n_slices)
        self.assert_(isfinite(log(
            vpif.gfn.compute_green_fn_from_distances(
                distances,
                denominator,
                start_slice, end_slice,
            )
        )))
    #@-node:gcross.20090918120753.2618:test_finite_log_half_path
    #@+node:gcross.20090918120753.2620:test_ignores_outside_distances
    @with_checker
    def test_ignores_outside_distances(self,
            n_slices = irange(6,102,4),
            denominator = unit_interval_float,
        ):
        start_slice = randint(1,n_slices-1)
        end_slice = randint(start_slice+1,n_slices)

        distances = rand(n_slices)

        distances[:start_slice-1] = 0
        distances[end_slice:start_slice] = 0

        self.assert_(
            vpif.gfn.compute_green_fn_from_distances(
                distances,
                denominator,
                start_slice, end_slice,
            ) > 0
        )
    #@-node:gcross.20090918120753.2620:test_ignores_outside_distances
    #@+node:gcross.20090918120753.2612:test_unit_interval
    @with_checker
    def test_unit_interval(self,
            n_slices = irange(6,102),
            magnitude = float,
            denominator = unit_interval_float,
        ):
        distances = rand(n_slices) * abs(magnitude)
        greens_fn = \
            vpif.gfn.compute_green_fn_from_distances(
                distances,
                denominator,
                1, n_slices,
            )
        self.assert_(greens_fn >= 0)
        self.assert_(greens_fn <= 1)
    #@-node:gcross.20090918120753.2612:test_unit_interval
    #@-others
#@-node:gcross.20090918120753.2604:compute_green_fn_from_distances
#@-others

tests = [
    gfn2_sp,
    gfn4_sp,
    gfn_hard_wall_contribution,
    compute_green_fn_from_distances,
    ]

if __name__ == "__main__":
    unittest.main()
#@nonl
#@-node:gcross.20090812093015.1742:@thin gfn.py
#@-leo
