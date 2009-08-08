#@+leo-ver=4-thin
#@+node:gcross.20090807144330.2253:@thin angular_momentum.py
#@@language python
#@@tabwidth -4

import unittest
from paycheck import with_checker
from paycheck.generator import positive_float, non_negative_float
from numpy import zeros, array, complex128
from numpy.linalg import norm
from tests import particle_paths_type
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
    def test_trivial_case(self):
        vector = array([2+0j],dtype=complex128)
        amplitudes = array([1+1j],dtype=complex128)
        vpi.angular_momentum.perform_special_matmul(vector,amplitudes)
        self.assertEqual(vector,array([2+2j],dtype=complex128))

    @with_checker([(complex,complex)])
    def test_random_cases(self,values):
        if(len(values) == 0): return
        vector, amplitudes = (array(values_,dtype=complex128) for values_ in zip(*values))
        original_vector = vector.copy()
        vpi.angular_momentum.perform_special_matmul(vector,amplitudes)
        partial_sum = 0
        for i in xrange(len(vector)):
            partial_sum += original_vector[i]
            self.assertEqual(vector[i],partial_sum*amplitudes[i])
#@-node:gcross.20090807144330.2254:perform_special_matmul
#@-others

tests = [
    get_rotation_plane_axes,
    perform_special_matmul
    ]

if __name__ == "__main__":
    unittest.main()
#@-node:gcross.20090807144330.2253:@thin angular_momentum.py
#@-leo
