#@+leo-ver=4-thin
#@+node:gcross.20090807144330.2151:@thin thermalize.py
#@@language python
#@@tabwidth -4

import unittest
from paycheck import with_checker
from paycheck.generator import positive_float, non_negative_float
from numpy import zeros, double
from numpy.linalg import norm
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
#@-others

tests = [
    accept_path
    ]

if __name__ == "__main__":
    unittest.main()
#@-node:gcross.20090807144330.2151:@thin thermalize.py
#@-leo
