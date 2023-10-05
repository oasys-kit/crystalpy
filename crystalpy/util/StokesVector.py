"""
Represents a Stokes vector (four components s0, s1, s2, s3). It accepts a stack (s0, etc. are arrays).
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
        self._s0 = numpy.array(element_list[0])
        self._s1 = numpy.array(element_list[1])
        self._s2 = numpy.array(element_list[2])
        self._s3 = numpy.array(element_list[3])

    def duplicate(self):
        """Duplicates a StokesVector.

        Returns
        -------
        StokesVector instance
            New StokesVector instance with identical x,y,z components.

        """
        return StokesVector(element_list=self.getList())


    def components(self):
        """Generates a numpy 1x4 array from the Stokes vector components.

        Returns
        -------
        numpy array
            The four stokes components.

        """
        return numpy.array(self.getList())

    @property
    def s0(self):
        """
        Gets s0 (first stokes parameter).

        Returns
        -------
        float or numpy array

        """
        return(self._s0)

    @property
    def s1(self):
        """
        Gets s1 (second stokes parameter).

        Returns
        -------
        float or numpy array

        """
        return(self._s1)

    @property
    def s2(self):
        """
        Gets s2 (first stokes parameter).

        Returns
        -------
        float or numpy array

        """
        return(self._s2)

    @property
    def s3(self):
        """
        Gets s3 (first stokes parameter).

        Returns
        -------
        float or numpy array

        """
        return(self._s3)

    def append(self, vector):
        """
        Appends a stoked-vector to the stack.

        Parameters
        ----------
        vector : instance of StokesVector

        """
        s00 = numpy.append(self.s0, vector.s0)
        s11 = numpy.append(self.s1, vector.s1)
        s22 = numpy.append(self.s2, vector.s2)
        s33 = numpy.append(self.s3, vector.s3)
        self.setFromArray([s00,s11,s22,s33])

    def concatenate(self, vector):
        """
        Concatenates a vector to the stack.

        Parameters
        ----------
        vector : instance of StokesVector

        Returns
        -------
        instance of StokesVector
            The resulting vector with the concatenation.

        """
        s00 = numpy.append(self.s0, vector.s0)
        s11 = numpy.append(self.s1, vector.s1)
        s22 = numpy.append(self.s2, vector.s2)
        s33 = numpy.append(self.s3, vector.s3)
        return StokesVector([s00,s11,s22,s33])

    def getS0(self):
        """Returns the S0 component.

        Returns
        -------

        float
            The S0 component.

        """
        return self.s0

    def getS1(self):
        """Returns the S1 component.

        Returns
        -------

        float
            The S1 component.

        """
        return self.s1

    def getS2(self):
        """Returns the S2 component.

        Returns
        -------

        float
            The S2 component.

        """
        return self.s2

    def getS3(self):
        """Returns the S3 component.

        Returns
        -------

        float
            The S3 component.

        """
        return self.s3

    def getList(self):
        """Generates a 1x4 list with the four Stokes components.

        Returns
        -------
        list
            list containing the Stokes parameters.

        """
        result = list()
        result.append(self._s0)
        result.append(self._s1)
        result.append(self._s2)
        result.append(self._s3)

        return result

    def setFromArray(self, array):
        """Set stokes components from a given array

        Parameters
        ----------
        array : list or numpy array


        """

        self._s0 = numpy.array(array[0])
        self._s1 = numpy.array(array[1])
        self._s2 = numpy.array(array[2])
        self._s3 = numpy.array(array[3])

    def setFromValues(self, s0, s1, s2, s3):
        """Set stokes components from given values

        Parameters
        ----------
        s0 : float
            
        s1 : float
            
        s2 : float
            
        s3 : float

        """

        self._s0 = numpy.array(s0)
        self._s1 = numpy.array(s1)
        self._s2 = numpy.array(s2)
        self._s3 = numpy.array(s3)

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
            return self._s3 / self._s0
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
        return "{S0} {S1} {S2} {S3}".format(S0=self._s0, S1=self._s1, S2=self._s2, S3=self._s3)

    def __eq__(self, candidate):
        if (self._s0 != candidate._s0).any():
            return False

        if (self._s1 != candidate._s1).any():
            return False

        if (self._s2 != candidate._s2).any():
            return False

        if (self._s3 != candidate._s3).any():
            return False

        return True

if __name__ == "__main__":
    element_list = [0.78177969457877930,
                         0.22595711869558588,
                         0.28797567756487550,
                         0.58551861060989900]
    stokes_vector = StokesVector(element_list)

    assert(stokes_vector.s0 == 0.78177969457877930)
    assert(stokes_vector.s1 == 0.22595711869558588)
    assert(stokes_vector.s2 == 0.28797567756487550)
    assert(stokes_vector.s3 == 0.58551861060989900)

    array1 = stokes_vector.components()
    array2 = stokes_vector.getList()

    assert( isinstance(array1, numpy.ndarray))
    assert( isinstance(array2, list))
    numpy.testing.assert_array_equal(array1, numpy.asarray(element_list))
    assert(array2 == element_list)

    pol_deg = stokes_vector.circularPolarizationDegree()
    assert( (isinstance(pol_deg, float) or (isinstance(pol_deg, numpy.ndarray))))
    assert( pol_deg ==  0.7489560226111716 )


    stokes_vector1 = StokesVector([[1, 1], [2, 2], [3, 3], [4, 4]])  # without final zeros
    stokes_vector2 = StokesVector([[1, 1], [2, 2], [3, 3], [4, 4]])  # without final zeros
    print(stokes_vector2.s0)
    assert(stokes_vector1 == stokes_vector1)
    assert(not(stokes_vector1 != stokes_vector2))
