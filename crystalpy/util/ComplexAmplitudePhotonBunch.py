"""
This object represents a bunch of "complex-amplitude"" photons. It contains a stack of ComplexAmplitudePhoton instances, characterized by photon energy, direction vector and sigma and pi complex amplitudes.
"""

import numpy as np

from crystalpy.util.PhotonBunch import PhotonBunch
from crystalpy.util.ComplexAmplitudePhoton import ComplexAmplitudePhoton
from crystalpy.util.Vector import Vector
import numpy

#todo: replace name "polarized" by "complex_amplitude"
class ComplexAmplitudePhotonBunchOld(PhotonBunch):
    """Constructor.

    Parameters
    ----------
    complex_amplitude_photons : list, optional
        List of ComplexAmplitudePhoton instances.

    """
    def __init__(self, complex_amplitude_photons=None):

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

#
# NEWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
#

class ComplexAmplitudePhotonBunch(ComplexAmplitudePhoton):
    """Constructor.

    New version that inheritates from ComplexAmplitudePhoton and uses stacks for more effcient stockage.

    Parameters
    ----------
    complex_amplitude_photons : list, optional
        List of ComplexAmplitudePhoton instances.

    """
    def __init__(self, complex_amplitude_photons=None):

        if complex_amplitude_photons == None:
            super().__init__(energy_in_ev=[], direction_vector=Vector([],[],[]), Esigma=[], Epi=[])
        else:
            n = len(complex_amplitude_photons)
            if n == 0:
                super().__init__(energy_in_ev=[], direction_vector=Vector([],[],[]), Esigma=[], Epi=[])
            else:
                energy = numpy.zeros(n)
                Esigma = numpy.zeros(n)
                Epi    = numpy.zeros(n)
                for el,i in enumerate(complex_amplitude_photons):
                    energy[i] = el.energy()
                    Esigma[i] = el.getComplexAmplitudeS()
                    Epi[i] = el.getComplexAmplitudeP()
                    vv = el.unitDirectionVector()
                    if i == 0:
                        v = Vector(
                            vv.components()[0],
                            vv.components()[1],
                            vv.components()[2],
                        )
                    else:
                        v = v.append(vv)
                self.setEnergy(energy)
                self.setUnitDirectionVector(v)
                super().__init__(energy_in_ev=energy, direction_vector=v, Esigma=Esigma, Epi=Epi)


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


    def addPhoton(self, to_be_added):
        """Adds a photon to the bunch.

        Parameters
        ----------
        to_be_added : Photon instance

        """
        self.setEnergy(numpy.append(self.energy(), to_be_added.energy()))
        self.setComplexAmplitudeS(numpy.append(self.getComplexAmplitudeS(), to_be_added.getComplexAmplitudeS()))
        self.setComplexAmplitudeP(numpy.append(self.getComplexAmplitudeP(), to_be_added.getComplexAmplitudeP()))
        self.setUnitDirectionVector(self.unitDirectionVector().concatenate(to_be_added.unitDirectionVector()))

    def addPhotonsFromList(self, to_be_added):
        """Adds a list of photons to the bunch.

        Parameters
        ----------
        to_be_added : list
            The photons to be added

        """
        for el in to_be_added:
            self.addPhoton(el)

    def addBunch(self, to_be_added):
        """Adds photons in a PhotonBunch instance.

        Parameters
        ----------
        to_be_added : PhotonBunch instance
            Photons to be added.


        """
        self.addPhoton(to_be_added)

    def getNumberOfPhotons(self):
        """Returns the number of photons in the bunch.

        Returns
        -------
        int
            Number of photons.
        """
        return self.unitDirectionVector().nStack()

    def getListOfPhotons(self):
        """Returns a list with the photons in the bunch.

        Returns
        -------
        list
            List with photons.
        """
        out = []
        v = self.unitDirectionVector()
        vx = v.components()[0]
        vy = v.components()[1]
        vz = v.components()[2]
        for i in range(self.getNumberOfPhotons()):
            out.append(ComplexAmplitudePhoton(
                              energy_in_ev=self.energy()[i],
                              direction_vector=Vector(vx[i], vy[i], vz[i]),
                              Esigma=self.getComplexAmplitudeS()[i],
                              Epi=self.getComplexAmplitudeP()[i],
                              )
                       )
        return out

    def getPhotonIndex(self, index):
        """Returns the photon in the bunch with a given index.

        Parameters
        ----------
        index : int
            The photon index to be referenced.

        Returns
        -------
        Photon instance
            The photon (referenced, not copied).

        """
        v = self.unitDirectionVector()
        vx = v.components()[0]
        vy = v.components()[1]
        vz = v.components()[2]
        return ComplexAmplitudePhoton(energy_in_ev=self.energy()[index],
                                      direction_vector=Vector(vx[index], vy[index], vz[index]),
                                      Esigma=self.getComplexAmplitudeS()[i],
                                      Epi=self.getComplexAmplitudeP()[i])

    def setPhotonIndex(self, index, polarized_photon):
        """Sets the photon in the bunch with a given index.

        Parameters
        ----------
        index : int
            The photon index to be modified.

        polarized_photon : Photon instance
            The photon to be stored.

        """
        energy = self.energy()
        v = self.unitDirectionVector()
        vx = v.components()[0]
        vy = v.components()[1]
        vz = v.components()[2]

        energy[index] = polarized_photon.energy()
        vx[index] = polarized_photon.unitDirectionVector().components()[0]
        vy[index] = polarized_photon.unitDirectionVector().components()[1]
        vz[index] = polarized_photon.unitDirectionVector().components()[2]
        self.setEnergy(energy)
        self.setUnitDirectionVector(Vector(vx, vy, vz))

    def keys(self):
        """return the keys of the dictionary resulting from toDictionary method"""
        return self.toDictionary().keys()

    def getArrayByKey(self, key):
        """Returns the array of a givem key in from toDictionary method

        Parameters
        ----------
        key :
            deviations', 's0', 's1', 's2', 's3'.

        Returns
        -------
        numpy array


        """
        return self.toDictionary()[key]

    def isMonochromatic(self, places):
        """Inquires about bunch monochromaticity.

        Parameters
        ----------
        places :
            number of decimal places to be taken into account for comparing energies.

        Returns
        -------
        bool
            True if all photons in the bunch have the same energy.

        """
        return self.energy().allEqual(self.energy()[0])

    def isUnidirectional(self):
        """Inquires if all photons in the bunch have the same direction.


        Returns
        -------
        bool
            True if all photons have the same direction.

        """
        first_direction = self.getPhotonIndex(0).unitDirectionVector()  # Vector object.

        # if the first element goes the same direction as all others, then all others share the same direction.
        for i in range(self.getNumberOfPhotons()):
            if first_direction != self.getPhotonIndex(i).unitDirectionVector():  # the precision is set to 7 decimal places.
                return False

        return True

    def __len__(self):
        return self.getNumberOfPhotons()

    def __iter__(self):
        return iter(self.getListOfPhotons())

    def energies(self):
        return self.energy()

    # def __getitem__(self, key):
    #     return self.polarized_photon_bunch[key]



if __name__ == "__main__":

    npoint = 10
    vx = numpy.zeros(npoint) + 0.0
    vy = numpy.zeros(npoint) + 1.0
    vz = numpy.zeros(npoint) + 0.0

    energy = numpy.zeros(npoint) + 3000.0

    photon_bunch1 = ComplexAmplitudePhotonBunch()
    photon_bunch2 = ComplexAmplitudePhotonBunch()

    photons_list = list()

    for i in range(npoint):
        photon = ComplexAmplitudePhoton(energy_in_ev=energy[i],
                        direction_vector=Vector(vx[i], vy[i], vz[i]),
                        Esigma=1+0j, Epi=0+1j)

        photon_bunch1.addPhoton(photon)
        photons_list.append(photon)

    print(photon_bunch1.toDictionary()['complexAmplitudeS'].shape)
    print(photon_bunch1.toDictionary()['complexAmplitudeP'].shape)
    print(photon_bunch1.toDictionary())