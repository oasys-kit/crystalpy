"""
Represents a photon defined by its energy, direction vector and pi and sigma complex amplitudes.
"""

from crystalpy.util.Photon import Photon
import numpy

class ComplexAmplitudePhoton(Photon):
    """Constructor.

    Parameters
    ----------
    energy_in_ev : float
        Photon energy in eV.

    direction_vector : Vector instance
        The direction of the photon (no need to be normalized).

    Esigma : complex
        The sigma-amplitude.

    Esigma : complex
        The pi-amplitude.

    """
    def __init__(self, energy_in_ev, direction_vector, Esigma=None, Epi=None):
        # Call base constructor.
        Photon.__init__(self, energy_in_ev, direction_vector)

        if Esigma is None:
            self._Esigma = numpy.array(1/numpy.sqrt(2)+0j)
        else:

            self._Esigma = numpy.array(Esigma, dtype=complex)

        if Epi is None:
            self._Epi = numpy.array(1/numpy.sqrt(2)+0j)
        else:
            self._Epi = numpy.array(Epi, dtype=complex)

    def rescaleEsigma(self, factor):
        """Multiply the sigma complex amplitude by a factor.

        Parameters
        ----------
        factor : float
            The multiplying factor.

        """
        self._Esigma *= numpy.array(factor, dtype=complex) # TODO: this cast may lose resolution


    def rescaleEpi(self, factor):
        """Multiply the pi complex amplitude by a factor.

        Parameters
        ----------
        factor : float
            The multiplying factor.

        """
        self._Epi *= numpy.array(factor, dtype=complex)  # TODO: this cast may lose resolution

    def getIntensityS(self):
        """Gets the sigma intensity.

        Returns
        -------
        float
            Intensity (sigma) of photon.

        """
        return numpy.abs(self._Esigma) ** 2

    def getIntensityP(self):
        """Gets the pi intensity.

        Returns
        -------
        float
            Intensity (pi) of photon.

        """
        return numpy.abs(self._Epi) ** 2

    def getIntensity(self):
        """Gets the total (sigma plus pi) intensity.

        Returns
        -------
        float
            Intensity of photon.

        """
        return self.getIntensityS() + self.getIntensityP()

    def getPhaseS(self):
        """Gets the sigma phase.

        Returns
        -------
        float
            Sigma-phase in radians.

        """
        return numpy.angle(numpy.array(self._Esigma, dtype=complex))

    def getPhaseP(self):
        """Gets the pi phase.

        Returns
        -------
        float
            Pi-phase in radians.

        """
        return numpy.angle(numpy.array(self._Epi, dtype=complex))

    def getComplexAmplitudeS(self):
        """Gets the sigma complex amplitude.

        Returns
        -------
        complex
            Sigma-complex amplitude.

        """
        return self._Esigma

    def getComplexAmplitudeP(self):
        """Gets the pi complex amplitude.

        Returns
        -------
        complex
            Pi-complex amplitude.

        """
        return self._Epi

    def setComplexAmplitudeS(self, value):
        """Sets the sigma complex amplitude.

        Parameters
        ----------
        value : complex, numpy array
            the value (complex) or numpy array (dtype=complex)

        """
        self._Esigma = value

    def setComplexAmplitudeP(self, value):
        """Sets the pi complex amplitude.

        Parameters
        ----------
        value : complex, numpy array
            the value (complex) or numpy array (dtype=complex)

        """
        self._Epi = value

    def __eq__(self, candidate):
        if ((self.energy() == candidate.energy() and
                self.unitDirectionVector() == candidate.unitDirectionVector()) and
                self._Esigma.complexAmplitude() == candidate._Esigma.complexAmplitude() and
                self._Epi.complexAmplitude() == candidate._Epi.complexAmplitude() ):
            return True

        return False

    def __ne__(self, candidate):
        return not (self == candidate)