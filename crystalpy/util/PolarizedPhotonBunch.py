"""
This object contains a list of PolarizedPhoton objects, characterized by energy, direction vector and Stokes vector.
This object is used as input to and output from the passive crystal widget.
"""
import numpy as np

from crystalpy.util.PolarizedPhoton import PolarizedPhoton
from crystalpy.util.PhotonBunch import PhotonBunch


# TODO: Create PhotonBunch with unpolarized photons
# TODO change methods to camelCase

class PolarizedPhotonBunch(PhotonBunch):
    """
    is a collection of PolarizedPhoton objects, making up the photon beam.
    """
    def __init__(self, polarized_photons=None):
        """
        :param polarized_photons: bunch of PolarizedPhoton objects.
        :type polarized_photons: list(PolarizedPhoton, PolarizedPhoton, ...)
        """
        if polarized_photons == None:
            self.polarized_photon_bunch = []
        else:
            self.polarized_photon_bunch = polarized_photons
        # self._set_dict()




    # # TODO: inheritated

    # def add(self, to_be_added):
    #     """
    #     :param to_be_added: PolarizedPhoton object(s) to be added to the bunch.
    #     :type to_be_added: PolarizedPhoton or list(PolarizedPhoton) of PhotonBunch
    #     """
    #     if type(to_be_added) == PolarizedPhoton:
    #         self.polarized_photon_bunch.append(to_be_added)
    #
    #     elif type(to_be_added) == list:
    #         self.polarized_photon_bunch.extend(to_be_added)
    #
    #     elif type(to_be_added) == PolarizedPhotonBunch:
    #         self.polarized_photon_bunch.extend(to_be_added.photon_bunch)
    #
    #     else:
    #         raise TypeError("The photon(s) could not be added to the bunch!")
    #
    #     self._set_dict()  # Update the array_dict attribute.


    # def __len__(self):
    #     return len(self.polarized_photon_bunch)
    #
    # def __iter__(self):
    #     return iter(self.polarized_photon_bunch)
    #
    # def __getitem__(self, key):
    #     return self.polarized_photon_bunch[key]
    # def get_number_of_photons(self):
    #     return len(self.polarized_photon_bunch)
    #
    # def get_list_of_photons(self):
    #     return self.polarized_photon_bunch
    #
    # def get_photon_index(self,index):
    #     return self.polarized_photon_bunch[index]
    #
    # def set_photon_index(self,index,polarized_photon):
    #
    #     self.polarized_photon_bunch[index] = polarized_photon
    # def is_monochromatic(self, places):
    #     """
    #     :param places: number of decimal places to be taken into account.
    #     :return: True if the bunch holds photons of the same energy.
    #     """
    #     first_energy = round(self.polarized_photon_bunch[0].energy(), places)
    #
    #     # if the first element has the same energy as all others, then all others share the same energy value.
    #     for polarized_photon in self:
    #         if first_energy != round(polarized_photon.energy(), places):
    #             return False
    #
    #     return True
    #
    # def is_unidirectional(self):
    #     """
    #     :return: True if the bunch holds photons going the same direction.
    #     """
    #     first_direction = self.polarized_photon_bunch[0].unitDirectionVector()  # Vector object.
    #
    #     # if the first element goes the same direction as all others, then all others share the same direction.
    #     for polarized_photon in self:
    #         if first_direction != polarized_photon.unitDirectionVector():  # the precision is set to 7 decimal places.
    #             return False
    #
    #     return True
    #
    # def get_array(self, key):
    #     """
    #     :param key: 'deviations', 's0', 's1', 's2', 's3'.
    #     :return: numpy.ndarray
    #     """
    #     return self.array_dict[key]


    def toDictionary(self):
        """
        defines a dictionary containing information about the bunch.
        """
        array_dict = PhotonBunch.toDictionary(self)

        stokes = np.zeros([4, len(self)])
        polarization_degrees = np.zeros(len(self))

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
        """
        :return: string containing the parameters characterizing each photon in the bunch.
        """
        bunch_string = str()
        for photon in self:
            string_to_attach = str(photon.energy()) + " " + \
                               photon.unitDirectionVector().toString() + "\n"
            bunch_string += string_to_attach
        return bunch_string


    # TODO create it "on the fly" do not store it in self
    # def _set_dict(self):
    #     """
    #     defines a dictionary containing information about the bunch.
    #     """
    #     self.array_dict = dict()
    #     energies = np.zeros(len(self))
    #     deviations = np.zeros(len(self))
    #     stokes = np.zeros([4, len(self)])
    #     directions = np.zeros([3, len(self)])
    #     polarization_degrees = np.zeros(len(self))
    #     i = 0
    #
    #     for polarized_photon in self:
    #         energies[i] = polarized_photon.energy()  # Photon.energy()
    #         deviations[i] = polarized_photon.deviation()
    #         stokes[0, i] = polarized_photon.stokesVector().s0
    #         stokes[1, i] = polarized_photon.stokesVector().s1
    #         stokes[2, i] = polarized_photon.stokesVector().s2
    #         stokes[3, i] = polarized_photon.stokesVector().s3
    #         directions[0, i] = polarized_photon.unitDirectionVector().components()[0]
    #         directions[1, i] = polarized_photon.unitDirectionVector().components()[1]
    #         directions[2, i] = polarized_photon.unitDirectionVector().components()[2]
    #         polarization_degrees[i] = polarized_photon.polarizationDegree()
    #         i += 1
    #
    #     self.array_dict["number of photons"] = i
    #     self.array_dict["energies"] = energies
    #     self.array_dict["deviations"] = deviations
    #     self.array_dict["s0"] = stokes[0, :]
    #     self.array_dict["s1"] = stokes[1, :]
    #     self.array_dict["s2"] = stokes[2, :]
    #     self.array_dict["s3"] = stokes[3, :]
    #     self.array_dict["vx"] = directions[0, :]
    #     self.array_dict["vy"] = directions[1, :]
    #     self.array_dict["vz"] = directions[2, :]
    #     self.array_dict["polarization degree"] = polarization_degrees
    #
    #
    #
    #
    #
    # def to_string(self):
    #     """
    #     :return: string containing the parameters characterizing each photon in the bunch.
    #     """
    #     bunch_string = str()
    #     for polarized_photon in self:
    #         string_to_attach = str(polarized_photon.energy()) + " " + \
    #                            polarized_photon.unitDirectionVector().toString() + " " + \
    #                            polarized_photon.stokesVector().toString() + " " + \
    #                            str(polarized_photon.polarizationDegree()) + "\n"
    #         bunch_string += string_to_attach
    #     return bunch_string
