"""
Unittest for Photon class.
"""

import unittest

from crystalpy.util.Photon import Photon
from crystalpy.util.Vector import Vector

import scipy.constants as codata
import numpy

class PhotonTest(unittest.TestCase):

    def testConstructorByDefault(self):
        photon = Photon()

        self.assertIsInstance(photon, Photon)
        self.assertEqual(photon.energy(), 1000.0)
        self.assertTrue(photon.unitDirectionVector() == Vector(0.0, 1.0, 0.0))

    def testConstructor(self):
        photon = Photon(4000, Vector(0, 0, 1))

        self.assertIsInstance(photon, Photon)
        self.assertEqual(photon.energy(), 4000)
        self.assertTrue(photon.unitDirectionVector() == Vector(0, 0, 1))

    def testEnergy(self):
        photon = Photon(4000, Vector(0, 0, 1))
        photon.setEnergy(8000.0)
        self.assertEqual(photon.energy(), 8000)

    def testWavelength(self):
        # Test data in eV : m.
        test_data = [   3,
                        4,
                        8,
                     5000,
                    10000, ]

        for energy in test_data:
            photon = Photon(energy, Vector(0, 0, 1))
            # print("Energy=%f, Wavelength=%f A (reference = %f A)"%(energy,1e10*photon.wavelength(),1e10*wavelength))
            wavelength = codata.h * codata.c / codata.e / energy
            self.assertAlmostEqual(photon.wavelength(),wavelength,places=1)

    def testWavenumber(self):
        # Test data in eV : m^-1.
        test_data = [   3,
                        4,
                        8,
                     5000,
                    10000,]

        for energy in test_data:
            photon = Photon(energy, Vector(0, 0, 1))
            wavelength = codata.h * codata.c / codata.e / energy
            wavenumber = 2 * numpy.pi / wavelength
            self.assertAlmostEqual(photon.wavenumber(),
                                   wavenumber, places=1)

    def testWavevector(self):
        direction = Vector(0, 0, 1)
        photon = Photon(5000.0, direction)

        wavevector = photon.wavevector()
        wavelength = codata.h * codata.c / codata.e / 5000
        aa = 2 * numpy.pi / wavelength
        self.assertAlmostEqual(wavevector.norm(),
                               aa, places=1)

        self.assertEqual(wavevector.getNormalizedVector(),
                         direction)

    def testUnitDirectionVector(self):
        photon = Photon(4000, Vector(0, 0, 5))

        self.assertTrue(photon.unitDirectionVector() == Vector(0, 0, 1))

    def testSetUnitDirectionVector(self):
        photon = Photon(4000, Vector(0, 0, 5))
        photon.setUnitDirectionVector(Vector(1,2,3))

        self.assertTrue(photon.unitDirectionVector() == Vector(1, 2, 3).getNormalizedVector())

    def testOperatorEqual(self):
        photon_one = Photon(4000, Vector(0, 0, 5))
        photon_two = Photon(4000, Vector(0, 1, 1))
        photon_three = Photon(2000, Vector(0, 0, 5))

        self.assertTrue(photon_one == photon_one)
        self.assertFalse(photon_one == photon_two)
        self.assertFalse(photon_one == photon_three)
        self.assertFalse(photon_two == photon_three)

    def testOperatorNotEqual(self):
        photon_one = Photon(4000, Vector(0, 0, 5))
        photon_two = Photon(4000, Vector(0, 1, 1))
        photon_three = Photon(2000, Vector(0, 0, 5))

        self.assertFalse(photon_one != photon_one)
        self.assertTrue(photon_one != photon_two)
        self.assertTrue(photon_one != photon_three)
        self.assertTrue(photon_two != photon_three)

    def testDuplicate(self):
        photon_one = Photon(4000, Vector(0, 0, 5))
        photon_two = photon_one.duplicate()


        self.assertTrue( photon_one == photon_two )

        photon_one.setEnergy(1000.0)
        self.assertFalse( photon_one == photon_two )