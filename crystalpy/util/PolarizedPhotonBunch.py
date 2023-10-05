"""
This object represents a bunch of polarized photons. It contains a stack of PolarizedPhoton instances, characterized by photon energy, direction vector and Stokes vector.
"""
import numpy

from crystalpy.util.PhotonBunch import PhotonBunch
from crystalpy.util.PolarizedPhoton import PolarizedPhoton
from crystalpy.util.Vector import Vector
from crystalpy.util.StokesVector import StokesVector


class PolarizedPhotonBunchOld(PhotonBunch):
    """The PolarizadPhotonBunch is is a collection of PolarizedPhoton objects, making up the polarized photon beam.

    Constructor.

    Parameters
    ----------
    polarized_photons : list, optional
        List of PolarizedPhoton instances.

    """


    def __init__(self, polarized_photons=None):
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




#
# NEWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
#
class PolarizedPhotonBunch(PolarizedPhoton):
    """The PolarizadPhotonBunch is is a PolarizedPhoton stack, making up the polarized photon beam.

    Constructor.

    Parameters
    ----------
    polarized_photons : list, optional
        List of PolarizedPhoton instances.

    """


    def __init__(self, polarized_photons=None):

        if polarized_photons == None:
            super().__init__(energy_in_ev=[],
                             direction_vector=Vector([],[],[]),
                             stokes_vector=StokesVector([[],[],[],[]]))
        else:
            n = len(polarized_photons)
            if n == 0:
                super().__init__(energy_in_ev=[],
                                 direction_vector=Vector([], [], []),
                                 stokes_vector=StokesVector([[], [], [], []]))
            else:
                energy = numpy.zeros(n)
                s0 = numpy.zeros(n)
                s1 = numpy.zeros(n)
                s2 = numpy.zeros(n)
                s3 = numpy.zeros(n)
                for el,i in enumerate(polarized_photons):
                    energy[i] = el.energy()
                    vv = el.unitDirectionVector()
                    ss = el.stokesVector()
                    if i == 0:
                        v = Vector(
                            vv.components()[0],
                            vv.components()[1],
                            vv.components()[2],
                        )
                        s = StokesVector([ss.s0, ss.s1, ss.s2, ss.s3])
                    else:
                        v.append(vv)
                        s.append(ss)
                self.setEnergy(energy)
                self.setUnitDirectionVector(v)
                super().__init__(energy_in_ev=energy, direction_vector=v, stokes_vector=s)



        # if polarized_photons == None:
        #     self.polarized_photon_bunch = []
        # else:
        #     self.polarized_photon_bunch = polarized_photons


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

    def addPhoton(self, to_be_added):
        """Adds a photon to the bunch.

        Parameters
        ----------
        to_be_added : Photon instance

        """
        self.setEnergy(numpy.append(self.energy(), to_be_added.energy()))
        self.setUnitDirectionVector(self.unitDirectionVector().concatenate(to_be_added.unitDirectionVector()))
        self.setStokesVector(self.stokesVector().concatenate(to_be_added.stokesVector()))


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
        s = self.stokesVector()
        for i in range(self.getNumberOfPhotons()):
            out.append(PolarizedPhoton(
                              energy_in_ev=self.energy()[i],
                              direction_vector=Vector(vx[i], vy[i], vz[i]),
                              stokes_vector=StokesVector([s.s0[i], s.s1[i], s.s2[i], s.s3[i]]),
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
        s = self.stokesVector()
        return PolarizedPhoton(energy_in_ev=self.energy()[index],
                                      direction_vector=Vector(vx[index], vy[index], vz[index]),
                                      stokes_vector=StokesVector([s.s0[index], s.s1[index], s.s2[index], s.s3[index]])
                               )

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
        s0 = self.stokesVector().s0
        s1 = self.stokesVector().s1
        s2 = self.stokesVector().s2
        s3 = self.stokesVector().s3

        energy[index] = polarized_photon.energy()
        vx[index] = polarized_photon.unitDirectionVector().components()[0]
        vy[index] = polarized_photon.unitDirectionVector().components()[1]
        vz[index] = polarized_photon.unitDirectionVector().components()[2]

        s0[index] = polarized_photon.stokesVector().s0
        s1[index] = polarized_photon.stokesVector().s1
        s2[index] = polarized_photon.stokesVector().s2
        s3[index] = polarized_photon.stokesVector().s3


        self.setEnergy(energy)
        self.setUnitDirectionVector(Vector(vx, vy, vz))
        self.setStokesVector(StokesVector([s0,s1,s2,s3]))

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