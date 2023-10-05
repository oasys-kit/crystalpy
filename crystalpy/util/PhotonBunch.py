"""
A stack of photons, each one characterized by energy and direction.
"""
import numpy
import copy
import scipy.constants as codata

from crystalpy.util.Vector import Vector
from crystalpy.util.Photon import Photon

#todo: delete
class PhotonBunchOld(object):
    """
    The PhotonBunch is Photon stack instances, making up the photon bunch or beam.

    Constructor.

    Parameters
    ----------
    photons : list
        List of Photon instances.

    """

    def __init__(self, photons=None):
        if photons == None:
            self.polarized_photon_bunch = []
        else:
            self.polarized_photon_bunch = photons

    @classmethod
    def initialize_from_energies_and_directions(cls, energies, V):
        """Construct a bunch from arrays with photon energies and directions

        Parameters
        ----------

        energies : list, numpy array
            
        V : Vector instance (with a tack of vectors)

        Returns
        -------
        PhotonBunch instance


        """
        if V.nStack() != energies.size:
            raise Exception("incompatible inputs")

        bunch = PhotonBunchOld()

        for i in range(energies.size):
            bunch.addPhoton(Photon(energy_in_ev=energies[i], direction_vector=V.extractStackItem(i)))

        return bunch

    def energies(self):
        """Return the energies of the photons.


        Returns
        -------
        numpy array
            The energies of the photons (copied, not referenced).

        """
        energies = numpy.zeros(len(self))
        for i,photon in enumerate(self):
            energies[i] = photon.energy()  # Photon.energy()
        return energies

    def energy(self): # just in case
        """Return the energies of the photons.


        Returns
        -------
        numpy array
            The energies of the photons (copied, not referenced).

        """
        return self.energies()

    def wavelength(self):
        """Return the wavelengths of the photons (in m).

        Returns
        -------
        numpy array
            The wavelengths of the photons.

        """
        E_in_Joule = self.energies() * codata.e # elementary_charge
        # Wavelength in meter
        wavelength = (codata.c * codata.h / E_in_Joule)
        return wavelength

    def wavenumber(self):
        """Return the wavenumbers of the photons (in m^-1).

        Returns
        -------
        numpy array
            The wavenumbers of the photons.

        """
        return (2.0 * numpy.pi) / self.wavelength()

    def unitDirectionVector(self):
        """Return the directions of the photons.

        Returns
        -------
        Vector instance
            The directions in stacked vectors.

        """
        X = numpy.zeros(len(self))
        Y = numpy.zeros(len(self))
        Z = numpy.zeros(len(self))
        for i,photon in enumerate(self):
            cc = photon.unitDirectionVector().components()
            X[i] = cc[0]
            Y[i] = cc[1]
            Z[i] = cc[2]
        return Vector.initializeFromComponents([X, Y, Z])


    def wavevector(self):
        """Return the wavevectors of the photons.

        Returns
        -------
        Vector instance
            The wavevectors in stacked vectors.

        """
        return self.unitDirectionVector().scalarMultiplication(self.wavenumber())

    def duplicate(self):
        """Return a clone of the PhotonBunch instance.

        Returns
        -------
        PhotonBunch instance

        """
        return copy.deepcopy(self)

    def setUnitDirectionVector(self, vector):
        """Sets the directions of the photons.

        Parameters
        ----------
        vector : Vector instance
            Stack of vectors with the directions.

        """
        for i,photon in enumerate(self):
            photon._unit_direction_vector = vector.extractStackItem(i)

    #
    # extend these methods when heritating from Photon
    #
    def toDictionary(self):
        """Created a dictionary containing information about the bunch.

        Returns
        -------
        dict
            Information in tags: "number of photons", "energies", "deviations", "vx", "vy" and "vz".

        """
        array_dict = dict()
        energies = numpy.zeros(len(self))
        deviations = numpy.zeros(len(self))
        directions = numpy.zeros([3, len(self)])

        i = 0

        for i,photon in enumerate(self):
            energies[i]      = photon.energy()  # Photon.energy()
            deviations[i]    = photon.deviation()
            directions[0, i] = photon.unitDirectionVector().components()[0]
            directions[1, i] = photon.unitDirectionVector().components()[1]
            directions[2, i] = photon.unitDirectionVector().components()[2]
            i += 1  # todo: very bizarre.... remove?

        array_dict["number of photons"] = i
        array_dict["energies"] = energies
        array_dict["deviations"] = deviations
        array_dict["vx"] = directions[0, :]
        array_dict["vy"] = directions[1, :]
        array_dict["vz"] = directions[2, :]

        return array_dict


    def toString(self):
        """Returns a string containing the parameters characterizing each photon in the bunch."""
        bunch_string = str()
        for photon in self:
            string_to_attach = str(photon.energy()) + " " + \
                               photon.unitDirectionVector().toString() + "\n"
            bunch_string += string_to_attach
        return bunch_string

    #
    # end of methods to be extended
    #

    def addPhoton(self, to_be_added):
        """Adds a photon to the bunch.

        Parameters
        ----------
        to_be_added : Photon instance

        """
        self.polarized_photon_bunch.append(to_be_added)


    def addPhotonsFromList(self, to_be_added):
        """Adds a list of photons to the bunch.

        Parameters
        ----------
        to_be_added : list
            The photons to be added

        """
        self.polarized_photon_bunch.extend(to_be_added)

    def addBunch(self, to_be_added):
        """Adds photons in a PhotonBunch instance.

        Parameters
        ----------
        to_be_added : PhotonBunch instance
            Photons to be added.
            

        """
        self.polarized_photon_bunch.extend(to_be_added.getListOfPhotons())

    def getNumberOfPhotons(self):
        """Returns the number of photons in the bunch.

        Returns
        -------
        int
            Number of photons.
        """
        return len(self.polarized_photon_bunch)

    def getListOfPhotons(self):
        """Returns a list with the photons in the bunch.

        Returns
        -------
        list
            List with photons.
        """
        return self.polarized_photon_bunch

    def getPhotonIndex(self,index):
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

        return self.polarized_photon_bunch[index]

    def setPhotonIndex(self,index,polarized_photon):
        """Sets the photon in the bunch with a given index.

        Parameters
        ----------
        index : int
            The photon index to be modified.

        polarized_photon : Photon instance
            The photon to be stored.

        """
        self.polarized_photon_bunch[index] = polarized_photon

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
        first_energy = round(self.polarized_photon_bunch[0].energy(), places)

        # if the first element has the same energy as all others, then all others share the same energy value.
        for photon in self:
            if first_energy != round(photon.energy(), places):
                return False

        return True

    def isUnidirectional(self):
        """Inquires if all photons in the bunch have the same direction.


        Returns
        -------
        bool
            True if all photons have the same direction.

        """
        first_direction = self.polarized_photon_bunch[0].unitDirectionVector()  # Vector object.

        # if the first element goes the same direction as all others, then all others share the same direction.
        for photon in self:
            if first_direction != photon.unitDirectionVector():  # the precision is set to 7 decimal places.
                return False

        return True

    def __len__(self):
        return len(self.polarized_photon_bunch)

    def __iter__(self):
        return iter(self.polarized_photon_bunch)

    def __getitem__(self, key):
        return self.polarized_photon_bunch[key]


#
#
#
class PhotonBunchDecorator(object):

    def energies(self):
        """Return the energies of the photons.


        Returns
        -------
        numpy array
            The energies of the photons (copied, not referenced).

        """
        return self.energy()


    def toString(self):
        """Returns a string table containing the energy and direction vector for each photon in the bunch."""
        bunch_string = str()
        for photon in self:
            string_to_attach = str(photon.energy()) + " (" + \
                               photon.unitDirectionVector().toString() + ")\n"
            bunch_string += string_to_attach
        return bunch_string

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
        for i in range(self.getNumberOfPhotons()):
            out.append(self.getPhotonIndex(i))
        return out

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

    #
    # these ones will be updated in ComplexAmplitudePhotonBunch and PolirizedPhotonBunch
    #
    def toDictionary(self):
        """Created a dictionary containing information about the bunch.

        Returns
        -------
        dict
            Information in tags: "number of photons", "energies", "deviations", "vx", "vy" and "vz".

        """
        array_dict = dict()
        e = self.energy()
        v = self.unitDirectionVector()
        n = v.nStack()

        array_dict["number of photons"] = n
        array_dict["energies"] = e
        array_dict["deviations"] = self.deviation()
        array_dict["vx"] = v.components()[0]
        array_dict["vy"] = v.components()[1]
        array_dict["vz"] = v.components()[2]

        return array_dict


    def addPhoton(self, to_be_added):
        """Adds a photon to the bunch.

        Parameters
        ----------
        to_be_added : Photon instance

        """
        self.setEnergy(numpy.append(self.energy(), to_be_added.energy()))
        self.setUnitDirectionVector(self.unitDirectionVector().concatenate(to_be_added.unitDirectionVector()))

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
        return Photon(energy_in_ev=self.energy()[index], direction_vector=Vector(vx[index], vy[index], vz[index]))

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

class PhotonBunch(Photon, PhotonBunchDecorator):
    """
    The PhotonBunch is a Photon stack.

    It inheritates from Photon and uses stacks for more efficient stockage. Additional methods
    useful for stacks or bunches are defined in PhotonBunchDecorator.

    Constructor.

    Parameters
    ----------
    photons : list
        List of Photon instances.

    """

    def __init__(self, photons=None):
        if photons == None:
            super().__init__(energy_in_ev=[], direction_vector=Vector([],[],[]))
        else:
            n = len(photons)
            energy = numpy.zeros(n)
            for el,i in enumerate(photons):
                energy[i] = el.energy()
                el.unitDirectionVector()
                if i == 0:
                    v = Vector(
                        el.components()[0],
                        el.components()[1],
                        el.components()[2],
                    )
                else:
                    v.append(el)
            self.setEnergy(energy)
            self.setUnitDirectionVector(v)
            super().__init__(energy_in_ev=energy, direction_vector=v)

    @classmethod
    def initialize_from_energies_and_directions(cls, energies, V):
        """Construct a bunch from arrays with photon energies and directions

        Parameters
        ----------

        energies : list, numpy array

        V : Vector instance (with a tack of vectors)

        Returns
        -------
        PhotonBunch instance


        """
        bunch = PhotonBunch()
        bunch.setEnergy(energies)
        bunch.setUnitDirectionVector(V)
        return bunch





if __name__ == "__main__":
    npoint = 10
    vx = numpy.zeros(npoint) + 0.0
    vy = numpy.zeros(npoint) + 1.0
    vz = numpy.zeros(npoint) + 0.0

    energy = numpy.zeros(npoint) + 3000.0

    photon_bunch1 = PhotonBunch()
    photon_bunch2 = PhotonBunch()

    photons_list = []

    for i in range(npoint):
        photon = Photon(energy_in_ev=energy[i],
                        direction_vector=Vector(vx[i], vy[i], vz[i]))

        photon_bunch1.addPhoton(photon)
        photons_list.append(photon)

    print(photon_bunch1.toDictionary())

    print(photon_bunch1.toString())