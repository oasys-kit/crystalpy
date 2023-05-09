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
    """
    is a collection of Photon objects, making up the photon beam.
    """
    def __init__(self, photons=None):
        """
        :param photons: bunch of Photon objects.
        :type photons: list(Photon, Photon, ...)
        """
        if photons == None:
            self.polarized_photon_bunch = []
        else:
            self.polarized_photon_bunch = photons
        # self._set_dict()

    @classmethod
    def initialize_from_energies_and_directions(cls, energies, V):
        if V.nStack() != energies.size:
            raise Exception("incompatible inputs")

        bunch = PhotonBunch()

        for i in range(energies.size):
            bunch.addPhoton(Photon(energy_in_ev=energies[i], direction_vector=V.extractStackItem(i)))

        return bunch

    def energies(self):
        energies = numpy.zeros(len(self))
        for i,photon in enumerate(self):
            energies[i]      = photon.energy()  # Photon.energy()
        return energies

    def energy(self): # just in case
        return self.energies()

    def wavelength(self):
        """
        :return: The photon wavelength in meter.
        """
        E_in_Joule = self.energies() * codata.e # elementary_charge
        # Wavelength in meter
        wavelength = (codata.c * codata.h / E_in_Joule)
        return wavelength

    def wavenumber(self):
        """
        :return: Wavenumber in m^-1.
        """
        return (2.0 * numpy.pi) / self.wavelength()

    def unitDirectionVector(self):
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
        """
        :return: Photon wavevector in m^-1.
        """
        return self.unitDirectionVector().scalarMultiplication(self.wavenumber())

    def duplicate(self):
        return copy.deepcopy(self)

    def setUnitDirectionVector(self, vector):
        for i,photon in enumerate(self):
            photon._unit_direction_vector = vector.extractStackItem(i)
            #
            # self._unit_direction_vector = vector.getNormalizedVector()

    #
    # extend these methods when heritating from Photon
    #
    def toDictionary(self):
        """
        defines a dictionary containing information about the bunch.
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
        """
        :return: string containing the parameters characterizing each photon in the bunch.
        """
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
        self.polarized_photon_bunch.append(to_be_added)


    def addPhotonsFromList(self, to_be_added):
        self.polarized_photon_bunch.extend(to_be_added)

    def addBunch(self, to_be_added):
        self.polarized_photon_bunch.extend(to_be_added.getListOfPhotons())

    def getNumberOfPhotons(self):
        return len(self.polarized_photon_bunch)

    def getListOfPhotons(self):
        return self.polarized_photon_bunch

    def getPhotonIndex(self,index):
        return self.polarized_photon_bunch[index]

    def setPhotonIndex(self,index,polarized_photon):
        self.polarized_photon_bunch[index] = polarized_photon

    def keys(self):
        return self.toDictionary().keys()

    def getArrayByKey(self, key):
        """
        :param key: 'deviations', 's0', 's1', 's2', 's3'.
        :return: numpy.ndarray
        """
        return self.toDictionary()[key]

    def isMonochromatic(self, places):
        """
        :param places: number of decimal places to be taken into account.
        :return: True if the bunch holds photons of the same energy.
        """
        first_energy = round(self.polarized_photon_bunch[0].energy(), places)

        # if the first element has the same energy as all others, then all others share the same energy value.
        for photon in self:
            if first_energy != round(photon.energy(), places):
                return False

        return True

    def isUnidirectional(self):
        """
        :return: True if the bunch holds photons going the same direction.
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