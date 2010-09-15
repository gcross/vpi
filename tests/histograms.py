#@+leo-ver=4-thin
#@+node:gcross.20090819083142.1375:@thin histograms.py
import unittest
from paycheck import *
from numpy import *
from numpy.random import rand
import __builtin__
import vpi.fortran as vpif

#@+others
#@+node:gcross.20090819083142.1376:class within_bins
class within_bins(unittest.TestCase):
    #@    @+others
    #@+node:gcross.20090819083142.1377:test_finite
    @with_checker
    def test_correct(self,bin=irange(-5,15),n_bins=irange(1,20)):
        if(bin >= 1 and bin <= n_bins):
            self.failUnless(vpif.histograms.within_bins(bin,n_bins))
        else:
            self.failIf(vpif.histograms.within_bins(bin,n_bins))
    #@nonl
    #@-node:gcross.20090819083142.1377:test_finite
    #@-others
#@-node:gcross.20090819083142.1376:class within_bins
#@+node:gcross.20090819093822.1380:class place_in_bin
class place_in_bin(unittest.TestCase):
    #@    @+others
    #@+node:gcross.20090819093822.1382:test_put_in_correct_bin_with_zero_offset
    @with_checker
    def test_put_in_correct_bin_with_zero_offset(self,n_bins=irange(1,10),number=unit_interval_float):
        histogram = zeros((n_bins),dtype='i',order='Fortran')
        vpif.histograms.place_in_bin(number,0,n_bins,histogram)
        correct_bin = floor(number * n_bins)
        for i in xrange(n_bins):
            if(i == correct_bin):
                self.assertEqual(1,histogram[i])
            else:
                self.assertEqual(0,histogram[i])
    #@-node:gcross.20090819093822.1382:test_put_in_correct_bin_with_zero_offset
    #@-others
#@-node:gcross.20090819093822.1380:class place_in_bin
#@-others

tests = [
    within_bins,
    place_in_bin,
    ]

if __name__ == "__main__":
    unittest.main()
#@nonl
#@-node:gcross.20090819083142.1375:@thin histograms.py
#@-leo
