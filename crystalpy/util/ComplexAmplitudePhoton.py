"""
This object contains a list of PolarizedPhoton objects, characterized by energy, direction vector and Stokes vector.
This object is used as input to and output from the passive crystal widget.
"""

from crystalpy.util.Photon import Photon
import numpy

# TODO create tests
class ComplexAmplitudePhoton(Photon):

    def __init__(self, energy_in_ev,direction_vector, Esigma=None,Epi=None):
        """
        Constructor.
        :param complex_amplitude: Complex amplitude of the wave.
        """

        # Call base constructor.
        Photon.__init__(self, energy_in_ev, direction_vector)

        # if Esigma == None:
        #     self._Esigma = ComplexAmplitude(1/numpy.sqrt(2)+0j)
        # else:
        #     if isinstance(Esigma,ComplexAmplitude):
        #         self._Esigma = Esigma
        #     else:
        #         self._Esigma = ComplexAmplitude(Esigma)
        #
        # if Epi == None:
        #     self._Epi = ComplexAmplitude(1/numpy.sqrt(2)+0j)
        # else:
        #     if isinstance(Epi,ComplexAmplitude):
        #         self._Epi = Epi
        #     else:
        #         self._Epi = ComplexAmplitude(Epi)

        if Esigma == None:
            self._Esigma = (1/numpy.sqrt(2)+0j)
        else:

            self._Esigma = Esigma

        if Epi == None:
            self._Epi = (1/numpy.sqrt(2)+0j)
        else:
            self._Epi = Epi



    def rescaleEsigma(self,factor):
        # if isinstance(factor,ComplexAmplitude):
        #     self._Esigma.rescale(factor.complexAmplitude())
        # else:
        #     self._Esigma.rescale(factor)

        self._Esigma *= factor


    def rescaleEpi(self,factor):
        # if isinstance(factor,ComplexAmplitude):
        #     self._Epi.rescale(factor.complexAmplitude())
        # else:
        #     self._Epi.rescale(factor)
        self._Epi *= factor

    def getIntensityS(self):
        """
        Sets the complex amplitude.
        :param complex_amplitude: Complex amplitude of the wave.
        """
        # return self._Esigma.intensity()
        return numpy.abs(self._Esigma) ** 2

    def getIntensityP(self):
        """
        Sets the complex amplitude.
        :param complex_amplitude: Complex amplitude of the wave.
        """
        # return self._Epi.intensity()
        return numpy.abs(self._Epi) ** 2

    def getIntensity(self):
        """
        Sets the complex amplitude.
        :param complex_amplitude: Complex amplitude of the wave.
        """
        return self.getIntensityS() + self.getIntensityP()

    def getPhaseS(self):
        return numpy.angle(numpy.array(self._Esigma, dtype=complex))

    def getPhaseP(self):
        return numpy.angle(numpy.array(self._Epi, dtype=complex))

    def getComplexAmplitudeS(self):
        return self._Esigma

    def getComplexAmplitudeP(self):
        return self._Epi

    # def duplicate(self):
    #     return ComplexAmplitudePhoton(self._energy_in_ev,
    #                            self._unit_direction_vector.duplicate(),
    #                            self._Esigma.complexAmplitude(),
    #                            self._Epi.complexAmplitude())

    def duplicate(self):
        return ComplexAmplitudePhoton(self._energy_in_ev,
                               self._unit_direction_vector.duplicate(),
                               self._Esigma,
                               self._Epi)

    def __eq__(self, candidate):
        """
        Determines if two polarized photons are identical (same energy, direction and polarization).
        :param candidate: Polarized photon to compare with.
        :return: True if equal otherwise False.
        """
        if ((self.energy() == candidate.energy() and
                self.unitDirectionVector() == candidate.unitDirectionVector()) and
                self._Esigma.complexAmplitude() == candidate._Esigma.complexAmplitude() and
                self._Epi.complexAmplitude() == candidate._Epi.complexAmplitude() ):
            return True

        return False

    # TODO not needed? inheritated?
    def __ne__(self, candidate):
        """
        Determines if two polarized photons are not identical (same energy, direction and polarization).
        :param candidate: Polarized photon to compare with.
        :return: True if not equal otherwise False.
        """
        return not (self == candidate)