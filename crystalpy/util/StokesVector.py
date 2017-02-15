"""
Represents a Stokes vector.
Except for energy all units are in SI. Energy is in eV.
"""
from numpy import asarray

# TODO change methods to camelCase

class StokesVector(object):

    #TODO define elements individually and not in list?

    def __init__(self, element_list=[0.0,0.0,0.0,0.0]):
        """
        Constructor.
        :param element_list: list containing the Stokes parameters S0,S1,S2,S3.
        """
        self.s0 = float(element_list[0])
        self.s1 = float(element_list[1])
        self.s2 = float(element_list[2])
        self.s3 = float(element_list[3])

    def duplicate(self):
        return StokesVector(self.get_array())

    # TODO change toList, toArray
    def get_array(self, numpy=True):
        """
        Generates a 1x4 array from the Stokes vector components.
        :param numpy: if True returns numpy.ndarray, if False returns list.
        :return: 1x4 array containing the Stokes parameters.
        """
        result = list()
        result.append(self.s0)
        result.append(self.s1)
        result.append(self.s2)
        result.append(self.s3)

        if numpy:
            return asarray(result)

        return result

    def set_array(self, array):

        self.s0 = float(array[0])
        self.s1 = float(array[1])
        self.s2 = float(array[2])
        self.s3 = float(array[3])

    def set_values(self, s0,s1,s2,s3):

        self.s0 = float(s0)
        self.s1 = float(s1)
        self.s2 = float(s2)
        self.s3 = float(s3)



    #TODO rename to circular_polarization_degree
    def polarization_degree(self):
        """
        Calculates the degree of circular polarization of the radiation
        described by the Stokes parameter.
        :return: degree of circular polarization
        """
        try:
            return self.s3 / self.s0
        except:
            return 0.0

    def to_string(self):
        """
        :return: a string object containing the four components of the Stokes vector.
        """
        return "{S0} {S1} {S2} {S3}".format(S0=self.s0, S1=self.s1, S2=self.s2, S3=self.s3)

    def __eq__(self, candidate):
        """
        Determines whether two Stokes vectors are equal.
        :param candidate: Stokes vector to compare to.
        :return: True if equal. False if not.
        """
        if self.s0 != candidate.s0:
            return False

        if self.s1 != candidate.s1:
            return False

        if self.s2 != candidate.s2:
            return False

        if self.s3 != candidate.s3:
            return False

        return True

    def __ne__(self, candidate):
        """
        Determines whether two Stokes vectors are not equal.
        :param candidate: Stokes vector to compare to.
        :return: True if not equal. False if equal.
        """
        return not self == candidate