"""
This object contains a list of PolarizedPhoton objects, characterized by energy, direction vector and Stokes vector.
This object is used as input to and output from the passive crystal widget.
"""
import numpy

from crystalpy.util.PhotonBunch import PhotonBunch


class PolarizedPhotonBunch(PhotonBunch):
    """The PolarizadPhotonBunch is is a collection of PolarizedPhoton objects, making up the polarized photon beam."""
    def __init__(self, polarized_photons=None):

        """Constructor.

        Parameters
        ----------
        polarized_photons : list
            List of PolarizedPhoton instances.

        Returns
        -------
        PolarizedPhotonBunch instance

        """

        if polarized_photons == None:
            self.polarized_photon_bunch = []
        else:
            self.polarized_photon_bunch = polarized_photons


    def toDictionary(self):
        """Created a dictionary containing information about the bunch.

        Returns
        -------
        dict
            Information in tags: "number of photons", "energies", "deviations", "vx", "vy", "vz", "s0", "s1", "s2", "s3" and "polarization degree".

        """

        """defines a dictionary containing information about the bunch."""
        array_dict = PhotonBunch.toDictionary(self)

        stokes = numpy.zeros([4, len(self)])
        polarization_degrees = numpy.zeros(len(self))

        for i,polarized_photon in enumerate(self):
            stokes[0, i] = polarized_photon.stokesVector().s0
            stokes[1, i] = polarized_photon.stokesVector().s1
            stokes[2, i] = polarized_photon.stokesVector().s2
            stokes[3, i] = polarized_photon.stokesVector().s3
            polarization_degrees[i] = polarized_photon.circularPolarizationDegree()

        array_dict["s0"] = stokes[0, :]
        array_dict["s1"] = stokes[1, :]
        array_dict["s2"] = stokes[2, :]
        array_dict["s3"] = stokes[3, :]
        array_dict["polarization degree"] = polarization_degrees

        return array_dict


    def toString(self):
        """Returns a string containing the parameters characterizing each polarized photon in the bunch."""
        bunch_string = str()

        for i in range(self.getNumberOfPhotons()):
            photon = self.getPhotonIndex(i)
            string_to_attach = str(photon.energy()) + " " + \
                               photon.unitDirectionVector().toString() + "\n"
            bunch_string += string_to_attach
        return bunch_string