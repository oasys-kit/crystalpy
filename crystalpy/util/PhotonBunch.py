"""
This object contains a list of Photon objects, characterized by energy and direction.

"""
import numpy
import scipy.constants as codata
from crystalpy.util.Vector import Vector
import copy
from crystalpy.util.Photon import Photon
from crystalpy.util.ComplexAmplitudePhoton import ComplexAmplitudePhoton
from crystalpy.util.PolarizedPhoton import PolarizedPhoton

class PhotonBunch(object):
    """The PhotonBunch is is a collection of Photon objects, making up the photon beam."""

    def __init__(self, photons=None):
        """Constructor.

        Parameters
        ----------
        photons : list
            List of Photon instances.

        Returns
        -------
        PhotonBunch instance

        """
        if photons == None:
            self.polarized_photon_bunch = []
        else:
            self.polarized_photon_bunch = photons

    @classmethod
    def initialize_from_energies_and_directions(cls, energies, V):
        """Construct a bunch from photon energies and directions


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

        bunch = PhotonBunch()

        for i in range(energies.size):
            bunch.addPhoton(Photon(energy_in_ev=energies[i], direction_vector=V.extractStackItem(i)))

        return bunch

    def energies(self):
        """Return the energies of the photons.


        Returns
        -------
        numpy array
            The energies of the photons.

        """
        energies = numpy.zeros(len(self))
        for i,photon in enumerate(self):
            energies[i]      = photon.energy()  # Photon.energy()
        return energies

    def energy(self): # just in case
        """Return the energies of the photons.


        Returns
        -------
        numpy array
            The energies of the photons.

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
            The referenced photon (not copied).

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