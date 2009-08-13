#@+leo-ver=4-thin
#@+node:gcross.20090812093015.1742:@thin gfn.py
#@@language python
#@@tabwidth -4

import unittest
from paycheck import with_checker, irange, frange
from numpy import zeros, double, dot
from numpy.linalg import norm
from numpy.random import rand
from random import randint
from tests import particle_paths_type
import vpi

#@+others
#@+node:gcross.20090812093015.1749:initialize_weights
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
#@+node:gcross.20090812093015.1743:gfn2_sp
class gfn2_sp(unittest.TestCase):

    @with_checker(irange(1,10),irange(1,5),float)
    def test_correctness(self,n_slices,n_particles,dt):
        sl_start = randint(1,n_slices)
        sl_end = randint(sl_start,n_slices)
        ip = randint(1,n_particles)
        U = rand(n_slices,n_particles)
        ln_gfn = vpi.gfn.gfn2_sp(sl_start,sl_end,ip,U,dt)
        self.assertAlmostEqual(ln_gfn,-dt*sum(U[sl_start-1:sl_end,ip-1]))
#@-node:gcross.20090812093015.1743:gfn2_sp
#@+node:gcross.20090812093015.1747:gfn4_sp
class gfn4_sp(unittest.TestCase):

    @with_checker(irange(1,15,2),irange(1,5),frange(0,1),frange(0,1))
    def test_correctness(self,c_slice,n_particles,lam,dt):
        n_slices = c_slice*2
        self.assertEqual(1,c_slice%2)

        (U_weight, gU2_weight) = initialize_weights(n_slices)

        sl_start = randint(1,n_slices)
        sl_end = randint(sl_start,n_slices)
        ip = randint(1,n_particles)
        U = rand(n_slices,n_particles)
        gradU2 = rand(n_slices)
        ln_gfn = vpi.gfn.gfn4_sp(sl_start,sl_end,ip,U,gradU2,U_weight,gU2_weight,lam,dt)

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
#@-node:gcross.20090812093015.1747:gfn4_sp
#@-others

tests = [
    gfn2_sp,
    gfn4_sp,
    ]

if __name__ == "__main__":
    unittest.main()
#@-node:gcross.20090812093015.1742:@thin gfn.py
#@-leo
