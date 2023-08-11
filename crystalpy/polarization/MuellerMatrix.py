"""
Represents a Mueller matrix.
See, e.g., https://en.wikipedia.org/wiki/Mueller_calculus
"""
import numpy

from crystalpy.util.StokesVector import StokesVector


class MuellerMatrix(object):
    """Constructor.

    Parameters
    ----------
    matrix :
        Matrix as a numpy array (4,4).

    """
    def __init__(self, matrix=numpy.zeros((4,4)) ):

        self.matrix = matrix

    @classmethod
    def initialize_as_general_linear_polarizer(cls,theta=0.0):
        """Creates a MuellerMatrix instance with a linear polarized.

        Parameters
        ----------
        theta : float, optional
            The the angle of the fast axis in rad. (Default value = 0.0)

        """
        mm = MuellerMatrix()
        mm.set_general_linear_polarizer(theta)
        return mm

    @classmethod
    def initialize_as_linear_polarizer_horizontal(cls):
        """Creates a MuellerMatrix instance with a horizontal linear polarized."""
        return cls.initialize_as_general_linear_polarizer(0.0)

    @classmethod
    def initialize_as_linear_polarizer_vertical(cls):
        """Creates a MuellerMatrix instance with a vertical linear polarized."""
        return cls.initialize_as_general_linear_polarizer(numpy.pi/2)

    @classmethod
    def initialize_as_linear_polarizer_plus45(cls):
        """Creates a MuellerMatrix instance with a +45 deg linear polarized."""
        return cls.initialize_as_general_linear_polarizer(numpy.pi/4)

    @classmethod
    def initialize_as_linear_polarizer_minus45(cls):
        """Creates a MuellerMatrix instance with a -45 deg linear polarized."""
        return cls.initialize_as_general_linear_polarizer(-numpy.pi/4)

    @classmethod
    def initialize_as_general_linear_retarder(cls,theta=0.0, delta=0.0):
        """Creates a MuellerMatrix instance with a phase retarder.

        Parameters
        ----------
        theta : float, optional
            The the angle of the fast axis in rad. (Default value = 0.0)
        delta : float, optional
            The phase difference between the fast and slow axis in rad. (Default value = 0.0)

        """
        mm = MuellerMatrix()
        mm.set_general_linear_retarder(theta,delta)
        return mm

    @classmethod
    def initialize_as_quarter_wave_plate_fast_vertical(cls):
        """Creare a MuellerMatrix instance with a quarter wave plate with fast axis vertical."""
        return cls.initialize_as_general_linear_retarder(numpy.pi/2,-numpy.pi/2)

    @classmethod
    def initialize_as_quarter_wave_plate_fast_horizontal(cls):
        """Creare a MuellerMatrix instance with a quarter wave plate with fast axis horizontal."""
        return cls.initialize_as_general_linear_retarder(0.0,-numpy.pi/2)

    @classmethod
    def initialize_as_half_wave_plate(cls):
        """Creare a MuellerMatrix instance with a half wave plate."""
        return cls.initialize_as_general_linear_retarder(0.0,numpy.pi)

    @classmethod
    def initialize_as_ideal_mirror(cls):
        """Creare a MuellerMatrix instance with a quarter wave plate with an ideal mirror."""
        return cls.initialize_as_general_linear_retarder(0.0,numpy.pi)

    @classmethod
    def initialize_as_filter(cls,transmission=1.0):
        """Creare a MuellerMatrix instance with a quarter wave plate with a filter.

        Parameters
        ----------
        transmission : float, optional
             The transmission value. (Default value = 1.0)

        Returns
        -------

        """
        return cls.initialize_as_general_linear_retarder(0.0,0.0).matrix_by_scalar(transmission)


    def from_matrix_to_elements(self, return_numpy=True):
        """Returns a list of flatten numpy array with the elements of a given matrix.
        If a list is needed one can use the numpy.array.tolist() method.

        Parameters
        ----------
        return_numpy : boolean, optional
            if True returns numpy.ndarray, if False returns list. (Default value = True)

        Returns
        -------
        numpy array or list
            [m00, m01, m02....mN0, mN1, mN2...]

        """
        matrix = numpy.asarray(self.matrix)
        result = matrix.flatten()

        if return_numpy:
            return result

        return list(result)

    def get_matrix(self):
        """Returns the muller matric (reference, not copied)."""
        return self.matrix

    def matrix_by_scalar(self, scalar):
        """Multiplies the matrix by a scalar.

        Parameters
        ----------
        scalar :
            the scalar factor.

        Returns
        -------
        MullerMatric instance
            the new Mueller matrix.

        """
        new_mueller_matrix = self.matrix * scalar

        return MuellerMatrix(new_mueller_matrix)

    def matrix_by_vector(self, vector, return_numpy=True):
        """Multiplies the matrix by a vector.

        Parameters
        ----------
        vector : numpy array
            the vector factor.
        return_numpy :
            if True returns numpy.ndarray, if False returns list. (Default value = True)

        Returns
        -------
        numpy array
            matrix * vector.

        """
        matrix = numpy.asarray(self.matrix)
        result = numpy.dot(matrix, vector)

        if return_numpy:
            return result

        return list(result)

    def vector_by_matrix(self, vector, return_numpy=True):
        """Multiplies a vector by the Muller matrix.

        Parameters
        ----------
        vector : numpy array
            the vector factor.
        return_numpy : boolean, optional
            if True returns numpy.ndarray, if False returns list. (Default value = True)

        Returns
        -------
        numpy array
            vector * matrix.

        """
        matrix = numpy.asarray(self.matrix)
        result = numpy.dot(vector, matrix)

        if return_numpy:
            return result

        return list(result)

    def mueller_times_mueller(self, matrix_2, mod=False):
        """Multiplies two Mueller matrices.

        Parameters
        ----------
        matrix_2 :
            Mueller matrix factor.
        mod : boolean, optional
            matrix multiplication is not commutative
            -> mod controls which of the two matrices is the first factor. (Default value = False, matrix_2 * mueller)

        Returns
        -------
        MuellerMatrix instance
            matrix * matrix_2 if mof=True
            matrix_2 * matrix if mof=false

        """
        matrix_1 = self.matrix

        if mod:
            product = numpy.dot(matrix_1, matrix_2)

        else:
            product = numpy.dot(matrix_2, matrix_1)

        return MuellerMatrix(product)

    def __eq__(self, candidate):
        for i in range(4):
            for j in range(4):

                if self.matrix[i, j] != candidate.matrix[i, j]:
                    return False

        return True

    def __ne__(self, candidate):
        return not self == candidate

    def set_general_linear_polarizer(self, theta):

        """Sets the Muller matrix as a linear polarizer. See [rf].

        Parameters
        ----------
        theta : float
            the angle of the fast axis in rad


        References
        ----------
        .. [rf] https://en.wikipedia.org/wiki/Mueller_calculus

        """

        # First row.
        self.matrix[0, 0] = 0.5
        self.matrix[0, 1] = 0.5 * numpy.cos(2*theta)
        self.matrix[0, 2] = 0.5 * numpy.sin(2*theta)
        self.matrix[0, 3] = 0.0

        # Second row.
        self.matrix[1, 0] = 0.5  * numpy.cos(2*theta)
        self.matrix[1, 1] = 0.5  * (numpy.cos(2*theta))**2
        self.matrix[1, 2] = 0.5  * numpy.sin(2*theta) * numpy.cos(2*theta)
        self.matrix[1, 3] = 0.0

        # Third row.
        self.matrix[2, 0] = 0.5 * numpy.sin(2*theta)
        self.matrix[2, 1] = 0.5 * numpy.sin(2*theta) * numpy.cos(2*theta)
        self.matrix[2, 2] = 0.5 * (numpy.sin(2*theta))**2
        self.matrix[2, 3] = 0.0

        # Fourth row.
        self.matrix[3, 0] = 0.0
        self.matrix[3, 1] = 0.0
        self.matrix[3, 2] = 0.0
        self.matrix[3, 3] = 0.0


    def set_general_linear_retarder(self, theta, delta=0.0):

        """Sets the Muller matrix as a generic line retarder. See [rf].

        Parameters
        ----------
        theta : float
            angle of fast axis in rad
        delta :
            phase difference in rad between the fast and slow axis in rad (Default value = 0.0)

        References
        ----------
        .. [rf] https://en.wikipedia.org/wiki/Mueller_calculus

        """

        # First row.
        self.matrix[0, 0] = 1.0
        self.matrix[0, 1] = 0.0
        self.matrix[0, 2] = 0.0
        self.matrix[0, 3] = 0.0

        # Second row.    (numpy.cos(2*theta))**2  (numpy.sin(2*theta))**2 (numpy.cos(delta))**2 (numpy.sin(delta))**2
        self.matrix[1, 0] = 0.0
        self.matrix[1, 1] = (numpy.cos(2*theta))**2 + numpy.cos(delta) * (numpy.sin(2*theta))**2
        self.matrix[1, 2] = numpy.cos(2*theta) * numpy.sin(2*theta) - numpy.cos(2*theta) * numpy.cos(delta) * numpy.sin(2*theta)
        self.matrix[1, 3] = numpy.sin(2*theta) * numpy.sin(delta)

        # Third row.
        self.matrix[2, 0] = 0.0
        self.matrix[2, 1] =  numpy.cos(2*theta) * numpy.sin(2*theta) - numpy.cos(2*theta) * numpy.cos(delta) * numpy.sin(2*theta)
        self.matrix[2, 2] =  numpy.cos(delta) * (numpy.cos(2*theta))**2 + (numpy.sin(2*theta))**2
        self.matrix[2, 3] = -numpy.cos(2*theta) * numpy.sin(delta)

        # Fourth row.
        self.matrix[3, 0] = 0.0
        self.matrix[3, 1] = -numpy.sin(2*theta) * numpy.sin(delta)
        self.matrix[3, 2] =  numpy.cos(2*theta) * numpy.sin(delta)
        self.matrix[3, 3] =  numpy.cos(delta)

    def calculate_stokes_vector(self, incoming_stokes_vector):
        """Takes an incoming Stokes vector, multiplies it by a Mueller matrix and gives an outgoing Stokes vector as a result.

        Parameters
        ----------
        incoming_stokes_vector : StokesVector instance
            The incoming vector.
            

        Returns
        -------
            StokesVector instance
                The resulting stokes vector.

        See Also
        --------
        crystalpy.util.StokesVector.StokesVector


        """
        # incoming_stokes_vector = self.incoming_stokes_vector.get_array()  # Stokes vector.
        element_list = self.matrix_by_vector(incoming_stokes_vector.getList())
        outgoing_stokes_vector = StokesVector(element_list)

        return outgoing_stokes_vector
