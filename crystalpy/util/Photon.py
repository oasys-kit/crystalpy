"""
Represents a photon defined by its energy and direction vector
"""

from crystalpy.util.Vector import Vector
import scipy.constants as codata
import numpy
import copy


class Photon(object):
    """Constructor.

    Parameters
    ----------
    energy_in_ev : float, optional
        Photon energy in eV.

    direction_vector : Vector instance, optional
        The direction of the photon (no need to be normalized).

    """
    def __init__(self, energy_in_ev=1000.0, direction_vector=Vector(0.0,1.0,0.0)):
        self._energy_in_ev = numpy.array(energy_in_ev)
        self._unit_direction_vector = direction_vector.getNormalizedVector()

        if self._unit_direction_vector.nStack() != self._energy_in_ev.size:
            raise Exception("Energy array must be of the same dimension of the vector stack.")

    def duplicate(self):
        """Return a clone of the Photon instance.

        Returns
        -------
        Photon instance

        """
        return copy.deepcopy(self)

    def energy(self):
        """Get the photon energy in eV.

        Returns
        -------
        float
            The photon energy.

        """
        return self._energy_in_ev

    def setEnergy(self,value):
        """

        Parameters
        ----------
        value : float
            The energy in eV to be set

        """
        self._energy_in_ev = value

    def wavelength(self):
        """Returns the photon wavelength in m.

        Returns
        -------
        float
            the photon wavelength

        """
        E_in_Joule = self.energy() * codata.e # elementary_charge
        # Wavelength in meter
        wavelength = (codata.c * codata.h / E_in_Joule)
        return wavelength

    def wavenumber(self):
        """Returns the photon wavenumber (2 pi / wavelength) in m^-1.

        Returns
        -------
        float
            the photon wavenumber

        """
        return (2.0 * numpy.pi) / self.wavelength()

    def wavevector(self):
        """Returns the photon wavevector in m^-1.

        Returns
        -------
        Vector instance
            the photon wavevector

        """
        return self.unitDirectionVector().scalarMultiplication(self.wavenumber())

    def unitDirectionVector(self):
        """Returns the photon direction vector.

        Returns
        -------
        Vector instance
            the photon dicection.

        """
        return self._unit_direction_vector

    def setUnitDirectionVector(self,vector=Vector(0,1,0)):
        """Sets the Photon direction.

        Parameters
        ----------
        vector : Vector instance
            The vector with the direction (may be not normalized).

        """
        self._unit_direction_vector = vector.getNormalizedVector()


    def deviation(self):
        """Returns the deviation angle (angle of the projection on YZ plane with Y axis)

        Returns
        -------
        float
            the deviation angle. The deviations are calculated supposing that the bunch moves along the y axis

        """
        vector = self.unitDirectionVector().components()  # ndarray([x, y, z])
        deviation = numpy.arctan2(vector[2], vector[1])

        return deviation

    def __eq__(self, candidate):
        if (self.energy() == candidate.energy() and
                self.unitDirectionVector() == candidate.unitDirectionVector()):
            return True

        return False

    def __ne__(self, candidate):
        return not (self == candidate)
