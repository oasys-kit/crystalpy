"""
This object contains a list of PolarizedPhoton objects, characterized by energy, direction vector and Stokes vector.
This object is used as input to and output from the passive crystal widget.
"""

from crystalpy.util.Photon import Photon
import numpy

# TODO create tests
class ComplexAmplitudePhoton(Photon):
    "This object represents a photon with energy, direction and complex amplitudes (sigma and pi)."
    def __init__(self, energy_in_ev, direction_vector, Esigma=None,Epi=None):
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

        Returns
        -------
            ComplexAmplitudePhoton instance.

        """

        # Call base constructor.
        Photon.__init__(self, energy_in_ev, direction_vector)

        if Esigma == None:
            self._Esigma = (1/numpy.sqrt(2)+0j)
        else:

            self._Esigma = Esigma

        if Epi == None:
            self._Epi = (1/numpy.sqrt(2)+0j)
        else:
            self._Epi = Epi



    def rescaleEsigma(self, factor):
        """Multiply the sigma complex amplitude by a factor.

        Parameters
        ----------
        factor : float
            The multiplying factor.

        """
        self._Esigma *= factor


    def rescaleEpi(self, factor):
        """Multiply the pi complex amplitude by a factor.

        Parameters
        ----------
        factor : float
            The multiplying factor.

        """
        self._Epi *= factor

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

    def duplicate(self):
        """Duplicates a complex-amplitude photon.

        Returns
        -------
        Photon instance
            New ComplexAmplitudePhoton instance with an identical photon.

        """
        return ComplexAmplitudePhoton(self._energy_in_ev,
                               self._unit_direction_vector.duplicate(),
                               self._Esigma,
                               self._Epi)

    def __eq__(self, candidate):
        if ((self.energy() == candidate.energy() and
                self.unitDirectionVector() == candidate.unitDirectionVector()) and
                self._Esigma.complexAmplitude() == candidate._Esigma.complexAmplitude() and
                self._Epi.complexAmplitude() == candidate._Epi.complexAmplitude() ):
            return True

        return False

    def __ne__(self, candidate):
        return not (self == candidate)