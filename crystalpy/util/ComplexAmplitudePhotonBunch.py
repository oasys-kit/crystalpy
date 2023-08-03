import numpy as np

from crystalpy.util.PhotonBunch import PhotonBunch


#todo: replace name "polarized" by "complex_amplitude"
class ComplexAmplitudePhotonBunch(PhotonBunch):
    """The ComplexAmplitudePhotonBunch is is a collection of ComplexAmplitudePhoton objects, making up the photon beam."""
    def __init__(self, complex_amplitude_photons=None):
        """Constructor.

        Parameters
        ----------
        complex_amplitude_photons : list
            List of ComplexAmplitudePhoton instances.

        Returns
        -------
        ComplexAmplitudePhotonBunch instance

        """
        if complex_amplitude_photons == None:
            self.polarized_photon_bunch = []
        else:
            self.polarized_photon_bunch = complex_amplitude_photons


    def rescaleEsigma(self, factor):
        """Multiply the sigma complex amplitudes by a factor array.

        Parameters
        ----------
        factor : list. numpy array
            The multiplying factor array.

        """
        for i, polarized_photon in enumerate(self):
            polarized_photon.rescaleEsigma(factor[i])

    def rescaleEpi(self, factor):
        """Multiply the pi complex amplitudes by a factor array.

        Parameters
        ----------
        factor : list. numpy array
            The multiplying factor array.

        """
        for i, polarized_photon in enumerate(self):
            polarized_photon.rescaleEpi(factor[i])

    def toDictionary(self):
        """Created a dictionary containing information about the bunch.

        Returns
        -------
        dict
            Information in tags: "number of photons", "energies", "deviations", "vx", "vy", "vz", "intensityS", "intensityP", "intensity", "phaseS", "phaseP", "complexAmplitudeS","complexAmplitudeP"

        """
        array_dict = PhotonBunch.toDictionary(self)

        intensityS = np.zeros(len(self))
        intensityP = np.zeros_like(intensityS)
        phaseS     = np.zeros_like(intensityS)
        phaseP     = np.zeros_like(intensityS)
        complexAmplitudeS = np.zeros_like(intensityS, dtype=complex)
        complexAmplitudeP = np.zeros_like(intensityS, dtype=complex)


        for i,polarized_photon in enumerate(self):
            intensityS[i] = polarized_photon.getIntensityS()
            intensityP[i] = polarized_photon.getIntensityP()
            phaseS    [i] = polarized_photon.getPhaseS()
            phaseP    [i] = polarized_photon.getPhaseP()
            complexAmplitudeS[i] = polarized_photon.getComplexAmplitudeS()
            complexAmplitudeP[i] = polarized_photon.getComplexAmplitudeP()


        array_dict["intensityS"] = intensityS
        array_dict["intensityP"] = intensityP
        array_dict["intensity"] = intensityS + intensityP
        array_dict["phaseS"] = phaseS
        array_dict["phaseP"] = phaseP
        array_dict["complexAmplitudeS"] = complexAmplitudeS
        array_dict["complexAmplitudeP"] = complexAmplitudeP


        return array_dict


    def toString(self):
        """Returns a string containing the parameters characterizing each photon in the bunch."""
        bunch_string = str()

        for i in range(self.getNumberOfPhotons()):
            photon = self.getPhotonIndex(i)
            string_to_attach = str(photon.energy()) + " " + \
                               photon.unitDirectionVector().toString() + "\n"
            bunch_string += string_to_attach
        return bunch_string