#@+leo-ver=4-thin
#@+node:gcross.20090914154930.2027:@thin erf.py
#@@language python
#@@tabwidth -4

import unittest
from paycheck import *
import scipy.special
import vpif

#@+others
#@+node:gcross.20090914154930.2028:erf
class erf(unittest.TestCase):

    @with_checker
    def test_correctness(self,x = float):
        self.assertAlmostEqual(scipy.special.erf(x),vpif.erfn.erf(x))
#@-node:gcross.20090914154930.2028:erf
#@+node:gcross.20090914154930.2035:erfc
class erfc(unittest.TestCase):

    @with_checker
    def test_correctness(self,x = float):
        self.assertAlmostEqual(scipy.special.erfc(x),vpif.erfn.erfc(x))
#@-node:gcross.20090914154930.2035:erfc
#@-others

tests = [
    erf,
    erfc,
    ]

if __name__ == "__main__":
    unittest.main()
#@nonl
#@-node:gcross.20090914154930.2027:@thin erf.py
#@-leo
