"""
Unittest for StokesVector class.
"""
import unittest

import numpy as np
from numpy.testing import assert_array_almost_equal

from crystalpy.util.StokesVector import StokesVector




class StokesVectorTest(unittest.TestCase):
    def setUp(self):
        self.element_list = [0.78177969457877930,
                             0.22595711869558588,
                             0.28797567756487550,
                             0.58551861060989900]
        self.stokes_vector = StokesVector(self.element_list)

    def testConstructor(self):
        self.assertIsInstance(self.stokes_vector, StokesVector)
        assert (self.stokes_vector.s0 == 0.78177969457877930)
        assert (self.stokes_vector.s1 == 0.22595711869558588)
        assert (self.stokes_vector.s2 == 0.28797567756487550)
        assert (self.stokes_vector.s3 == 0.58551861060989900)

    def testGetArray(self):
        array1 = self.stokes_vector.components()
        array2 = self.stokes_vector.getList()

        assert (isinstance(array1, np.ndarray))
        assert (isinstance(array2, list))
        np.testing.assert_array_equal(array1, np.asarray(self.element_list))
        assert (array2 == self.element_list)

    def testPolarizationDegree(self):
        pol_deg = self.stokes_vector.circularPolarizationDegree()
        assert ((isinstance(pol_deg, float) or (isinstance(pol_deg, np.ndarray))))
        assert (pol_deg == 0.7489560226111716)

    def test_operator_equal(self):
        stokes_vector1 = self.stokes_vector
        stokes_vector2 = StokesVector([0.7817796945787793,
                                       0.22595711869558588,
                                       0.2879756775648755,
                                       0.585518610609899])  # without final zeros
        stokes_vector3 = StokesVector([round(0.78177969457877930, 6),
                                       round(0.22595711869558588, 6),
                                       round(0.28797567756487550, 6),
                                       round(0.58551861060989900, 6)])  # rounded float
        self.assertTrue(stokes_vector1 == stokes_vector1)  # identity
        self.assertTrue(stokes_vector1 == stokes_vector2)
        self.assertFalse(stokes_vector1 == stokes_vector3)

    def test_operator_not_equal(self):
        stokes_vector1 = self.stokes_vector
        stokes_vector2 = StokesVector([0.7817796945787793,
                                       0.22595711869558588,
                                       0.2879756775648755,
                                       0.585518610609899])  # without final zeros
        stokes_vector3 = StokesVector([round(0.78177969457877930, 6),
                                       round(0.22595711869558588, 6),
                                       round(0.28797567756487550, 6),
                                       round(0.58551861060989900, 6)])  # rounded float
        self.assertFalse(stokes_vector1 != stokes_vector1)
        self.assertFalse(stokes_vector1 != stokes_vector2)
        self.assertTrue(stokes_vector1 != stokes_vector3)

    def testDuplicate(self):
        v1 = StokesVector([1,2,3,4])
        v2 = v1.duplicate()

        self.assertTrue(v1 == v2)

        v1.setFromValues(0,2,3,4)

        self.assertFalse(v1 == v2)

    def testArray(self):
        stokes_vector1 = StokesVector([[1, 1], [2, 2], [3, 3], [4, 4]])  # without final zeros
        stokes_vector2 = StokesVector([[1, 1], [2, 2], [3, 3], [4, 4]])  # without final zeros

        assert (stokes_vector1 == stokes_vector1)
        assert (not (stokes_vector1 != stokes_vector2))
