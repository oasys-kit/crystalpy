"""
Represents a Stokes vector.

Except for energy all units are in SI. Energy is in eV.
"""
import numpy

class StokesVector(object):
    """StokesVector Constructor.

    Parameters
    ----------
    element_list : list, optional
        the Stokes parameters [S0,S1,S2,S3]

    """
    def __init__(self, element_list=[0.0,0.0,0.0,0.0]):
        self.s0 = float(element_list[0])
        self.s1 = float(element_list[1])
        self.s2 = float(element_list[2])
        self.s3 = float(element_list[3])

    def duplicate(self):
        """Duplicates a StokesVector.

        Returns
        -------
        StokesVector instance
            New StokesVector instance with identical x,y,z components.

        """
        return StokesVector(self.components())


    def components(self):
        """Generates a numpy 1x4 array from the Stokes vector components.

        Returns
        -------
        numpy array
            The four stokes components.

        """
        return numpy.array(self.getList())

    def getS0(self):
        """Returns the S0 component.

        Returns
        -------

        float
            The S0 component.

        """
        return self.components()[0]

    def getS1(self):
        """Returns the S1 component.

        Returns
        -------

        float
            The S1 component.

        """
        return self.components()[1]

    def getS2(self):
        """Returns the S2 component.

        Returns
        -------

        float
            The S2 component.

        """
        return self.components()[2]

    def getS3(self):
        """Returns the S3 component.

        Returns
        -------

        float
            The S3 component.

        """
        return self.components()[3]

    def getList(self):
        """Generates a 1x4 list with the four Stokes components.

        Returns
        -------
        list
            list containing the Stokes parameters.

        """
        result = list()
        result.append(self.s0)
        result.append(self.s1)
        result.append(self.s2)
        result.append(self.s3)

        return result

    def setFromArray(self, array):
        """Set stokes components from a given array

        Parameters
        ----------
        array : list or numpy array


        """

        self.s0 = float(array[0])
        self.s1 = float(array[1])
        self.s2 = float(array[2])
        self.s3 = float(array[3])

    def setFromValues(self, s0,s1,s2,s3):
        """Set stokes components from given values

        Parameters
        ----------
        s0 : float
            
        s1 : float
            
        s2 : float
            
        s3 : float

        """

        self.s0 = float(s0)
        self.s1 = float(s1)
        self.s2 = float(s2)
        self.s3 = float(s3)

    def circularPolarizationDegree(self):
        """Calculates the degree of circular polarization of the radiation described by the Stokes parameter.

        Parameters
        ----------

        Returns
        -------
        float
            Degree of circular polarization S3/S0

        """
        try:
            return self.s3 / self.s0
        except:
            return 0.0

    def toString(self):
        """Returns a string with the four Stokes parameters (separated by a blank).


        Returns
        -------
        str
            the four Stokes parameters.

        """

        """:return: a string object containing the four components of the Stokes vector."""
        return "{S0} {S1} {S2} {S3}".format(S0=self.s0, S1=self.s1, S2=self.s2, S3=self.s3)

    def __eq__(self, candidate):
        if self.s0 != candidate.s0:
            return False

        if self.s1 != candidate.s1:
            return False

        if self.s2 != candidate.s2:
            return False

        if self.s3 != candidate.s3:
            return False

        return True

