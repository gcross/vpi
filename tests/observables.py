#@+leo-ver=4-thin
#@+node:gcross.20090818081913.1356:@thin observables.py

import unittest
from paycheck import *
from numpy import *
from numpy.random import rand
import __builtin__
import vpif

#@+others
#@+node:gcross.20090818081913.1357:class compute_local_energy_estimate
class compute_local_energy_estimate(unittest.TestCase):
    #@    @+others
    #@+node:gcross.20090818081913.1358:test_finite
    @with_checker
    def test_finite(self,n_particles=irange(1,10),n_dimensions=irange(1,5),lambda_=unit_interval_float):
        U = rand(n_particles)
        gradient_of_log_trial_fn = array(rand(n_particles,n_dimensions),dtype=double,order='Fortran')
        laplacian_of_log_trial_fn = rand(1)
        energy = vpif.observables.compute_local_energy_estimate(
            U,
            gradient_of_log_trial_fn, laplacian_of_log_trial_fn,
            lambda_
            )
        self.assert_(isfinite(energy))
    #@nonl
    #@-node:gcross.20090818081913.1358:test_finite
    #@-others
#@-node:gcross.20090818081913.1357:class compute_local_energy_estimate
#@-others

tests = [
    compute_local_energy_estimate,
    ]

if __name__ == "__main__":
    unittest.main()
#@nonl
#@-node:gcross.20090818081913.1356:@thin observables.py
#@-leo
