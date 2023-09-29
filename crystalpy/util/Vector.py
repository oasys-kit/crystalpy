"""
Represents a 3d vector. Accept arrays of vectors.

vector.components() gives a numpy array with shape (3) or (3, npoints)
"""
import numpy


class Vector(object):
    """Vector Constructor.

    Parameters
    ----------
    x : float
        x component. It can be an array.
    y : float
        y component. It can be an array of the same size of x.
    z : float
        z component. It can be an array of the same size of x.

    """
    def __init__(self, x, y, z):
        self.setComponents(x, y, z)

    @staticmethod
    def initializeFromComponents(components):
        """Creates a vector from a list/array of at least three elements.

        Parameters
        ----------
        components : list, numpy array
            A list with the [x,y,z] components of the vector.


        """
        return Vector(components[0],
                      components[1],
                      components[2])

    def duplicate(self):
        """Duplicates a vector.

        Returns
        -------
        result : Vector instance
            New Vector instance with identical x,y,z components.

        """
        return Vector(self.components()[0],self.components()[1],self.components()[2])

    def setComponents(self, x, y, z):
        """Sets vector components.

        Parameters
        ----------
        x :
            x component.
        y :
            y component.
        z :
            z component.

        """
        self._components = numpy.asarray([x, y, z])

    def components(self):
        """Returns the components of this vector a 1d three element array.

        Returns
        -------
        numpy array
            The Vector components (referenced, not copied).

        """
        return self._components

    def componentsStack(self):
        """Returns the components stack of shape (3, npoints) of this vector.

        Returns
        -------
        numpy array
            The Vector components (referenced, not copied).

        """
        if self.isArray():
            return self.components()
        else:
            return self.components().reshape((-1,1))


    def nStack(self):
        """Returns the number of stacked vectors

        Returns
        -------
        int
            Number of stacked vectors

        """
        s = numpy.array(self.components().shape)
        if s.size == 1:
            return 1
        else:
            return s[1]

    def extractStackItem(self, i):
        """

        Parameters
        ----------
        i : int
            Index of vector to be extracted
            

        Returns
        -------
        Vector instance
            A new vector instance with components of the i-th element only.

        """
        x = self.getX()
        y = self.getY()
        z = self.getZ()
        return Vector.initializeFromComponents([x[i], y[i], z[i]])

    def isArray(self):
        """Inform if the Vector instance contains a stack of vectors.

        Returns
        -------
        bool
            True if Vector contsins a stack of vector, False if there is a single vector.

        """
        if self.nStack() == 1:
            return False
        else:
            return True

    def getX(self):
        """Returns the x component.

        Returns
        -------
        float, numpy array
            The x component(s) (referenced, not copied).

        """
        return self.components()[0]

    def getY(self):
        """Returns the y component (referenced, not copied).

        Returns
        -------

        float, numpy array
            The y component(s).

        """
        return self.components()[1]

    def getZ(self):
        """Returns the z component (referenced, not copied).

        Returns
        -------

        float, numpy array
            The z component(s).

        """
        return self.components()[2]

    def addVector(self, summand):
        """Adds two vectors.

        Parameters
        ----------
        summand : Vector instance
            The vector to be added.

        Returns
        -------
        Vector instance
            The sum as a new vector.

        """

        wX = self.getX() + summand.getX()
        wY = self.getY() + summand.getY()
        wZ = self.getZ() + summand.getZ()
        return Vector.initializeFromComponents([wX, wY, wZ])


    def scalarMultiplication(self, k):
        """Scalar multiplication of a vectors.

        Parameters
        ----------
        k :
            The scalar to multiply with.

        Returns
        -------
        Vector instance
            The scalar multiplied vector as a new vector.

        """
        return Vector.initializeFromComponents([self.getX() * k, self.getY() * k, self.getZ() * k])


    def subtractVector(self, tosubstract):
        """Subtract a vector from this Vector instance.

        Parameters
        ----------

        tosubstract :
            Vector to subtract.

        Returns
        -------
        Vector instance
            The difference of the two vectors as a new vector.

        """
        result = self.addVector(tosubstract.scalarMultiplication(-1.0))
        return result

    def scalarProduct(self, factor):
        """Calculates the scalar product of this vector with the given vector.

        Parameters
        ----------
        factor :
            The vector to calculate the scalar product with.

        Returns
        -------
        Vector instance
            Scalar product of the two vectors as a new vector.

        """
        # scalar_product = numpy.dot(self.components(), factor.components())
        # scalar_product = numpy.sum( self.components() * factor.components(), axis=0)
        # return scalar_product

        wX = self.getX() * factor.getX()
        wY = self.getY() * factor.getY()
        wZ = self.getZ() * factor.getZ()
        return wX + wY + wZ


    def crossProduct(self, factor):
        """Calculates the cross product of two vectors.

        Parameters
        ----------
        factor :
            The vector to form the cross product with.

        Returns
        -------
        Vector instance
            Cross product of the two vectors as a new vector.

        """
        uX = self.getX()
        uY = self.getY()
        uZ = self.getZ()
        vX = factor.getX()
        vY = factor.getY()
        vZ = factor.getZ()

        wX = uY * vZ - uZ * vY
        wY = uZ * vX - uX * vZ
        wZ = uX * vY - uY * vX

        return Vector.initializeFromComponents([wX, wY, wZ])


    def norm(self):
        """Returns the standard norm of this norm.
        :return: Norm of this vector,

        Returns
        -------
        float, numpy array
            The vector norm.

        """
        norm = self.scalarProduct(self) ** 0.5
        return norm

    def getNormalizedVector(self):
        """Returns a normalized vector of this vector.
        :return: Normalized vector of this vector.

        Returns
        -------
        Vector instance
            The Normalized vector as a new vector.

        """
        return self.scalarMultiplication(self.norm() ** -1.0)

    def rotateAroundAxis(self, rotation_axis, angle):
        """Rotates the vector around an axis. It uses the Rodrigues formula [rf]_

        Parameters
        ----------
        rotation_axis : Vector instance
            Vector specifying the rotation axis (not necessarily unit vector).

        angle : float
            Rotation angle in radiants.

        Returns
        -------
        Vector instance
            Rotated vector as a new vector.

        References
        ----------
        .. [rf] http://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula

        """
        unit_rotation_axis = rotation_axis.getNormalizedVector()

        rotated_vector = self.scalarMultiplication(numpy.cos(angle))

        tmp_vector = unit_rotation_axis.crossProduct(self)
        tmp_vector = tmp_vector.scalarMultiplication(numpy.sin(angle))
        rotated_vector = rotated_vector.addVector(tmp_vector)

        scalar_factor = self.scalarProduct(unit_rotation_axis) * (1.0 - numpy.cos(angle))
        tmp_vector = unit_rotation_axis.scalarMultiplication(scalar_factor)
        rotated_vector = rotated_vector.addVector(tmp_vector)

        return rotated_vector

    def parallelTo(self, vector):
        """Returns the (parallel) projection of this vector along the given vector.

        Parameters
        ----------
        vector : Vector instance
            Vector defining the projection direction. It does not need to be normalized.

        Returns
        -------
        Vector instance
            Parallel projection along the vector.

        """
        unit_direction = vector.getNormalizedVector()
        projection_in_direction = self.scalarProduct(unit_direction)
        parallel_projection = unit_direction.scalarMultiplication(projection_in_direction)

        return parallel_projection

    def perpendicularTo(self, vector):
        """Returns the projection perpendicular to the given vector.

        Parameters
        ----------
        vector : Vector instance
            Vector that defines the direction. It does not need to be normalized.

        Returns
        -------
        Vector instance
            Projection perpendicular to the given vector.

        """
        perpendicular = self.subtractVector(self.parallelTo(vector))
        return perpendicular

    def getOnePerpendicularVector(self):
        """Returns one arbitrary vector perpendicular to this vector.

        Returns
        -------
        Vector instance
            One arbitrary vector perpendicular to this vector.

        """
        n = self.nStack()
        if n == 1:
            vector_y = Vector(0, 1, 0)
            vector_z = Vector(0, 0, 1)
        else:
            vector_y = Vector.initializeFromComponents( [numpy.zeros(n), numpy.ones(n), numpy.zeros(n)])
            vector_z = Vector.initializeFromComponents( [numpy.zeros(n), numpy.zeros(n), numpy.ones(n)])

        if self.getNormalizedVector() == vector_z:
            return vector_y

        vector_perpendicular = vector_z.perpendicularTo(self)
        vector_perpendicular = vector_perpendicular.getNormalizedVector()

        return vector_perpendicular

    def angle(self, factor):
        """Return the angle between this vector and the given vector.

        Parameters
        ----------
        factor : Vector instance
            Vector to determine the angle with.

        Returns
        -------
        float or numpy array
            Angle between this vector and the given vector.

        """
        n1 = self.getNormalizedVector()
        n2 = factor.getNormalizedVector()

        # Determine angle between the two vectors.
        cos_angle = n1.scalarProduct(n2)
        cos_angle= numpy.clip(cos_angle, -1, 1) # just in case...
        angle = numpy.arccos(cos_angle)

        return angle

    def scattering_on_surface(self, NORMAL, H):
        """Returns K_OUT vector following the scattering equation at a surface:
            K_OUT_parallel = K_IN_parallel + H_parallel
            |K_OUT| = |K_IN|

        Parameters
        ----------
        NORMAL : instance of Vector
            The vector normal to the surface.
        H : instance of Vector
            The scattering vector.

        Returns
        -------
        Vector instance
            Vector with K_OUT.

        """
        H_perp = NORMAL.scalarMultiplication(H.scalarProduct(NORMAL))
        H_par = H.subtractVector(H_perp)

        K_IN_perp = NORMAL.scalarMultiplication(self.scalarProduct(NORMAL))
        K_IN_par = self.subtractVector(K_IN_perp)

        K_OUT_par = K_IN_par.addVector(H_par)
        K_OUT_perp = NORMAL.scalarMultiplication(
            numpy.sqrt(self.norm() ** 2 - K_OUT_par.norm() ** 2))
        K_OUT = K_OUT_par.addVector(K_OUT_perp)
        return K_OUT

    def getVectorWithAngle(self, angle):
        """Returns one arbitrary vector with the vector rotated with a given angle.

        Parameters
        ----------
        angle : float
            The requested angle in radiants.

        Returns
        -------
        Vector instance
            Vector with given angle to this vector.

        """
        vector_perpendicular = self.getOnePerpendicularVector()
        vector_with_angle = self.rotateAroundAxis(vector_perpendicular, angle)

        return vector_with_angle

    def printComponents(self):
        """Print vector components
        """
        print(self.toString())

    def toString(self):
        """Returns a string object containing the components of the vector.

        Returns
        -------
        str
            Components of this vector.

        """
        return "{Vx} {Vy} {Vz}".format(Vx=self.components()[0],
                                       Vy=self.components()[1],
                                       Vz=self.components()[2])


    def __eq__(self, candidate):
        return numpy.linalg.norm(self.components()
                              -
                              candidate.components()) < 1.e-7

    def __ne__(self, candidate):
        return not (self == candidate)

    def __add__(self, o):
        return self.addVector(o)

    def __sub__(self, o):
        return self.subtractVector(o)

    def __mul__(self, o):
        return self.scalarMultiplication(o)


if __name__ == "__main__":
    vector = Vector(1, 2, 3)
    print(vector.components())
    print("dot: ", (vector.scalarProduct(vector)))
    print("shape: ", vector.components().shape)
    print("v x v =", vector.crossProduct(vector).components())
    print("3 v : ", vector.scalarMultiplication(3).components(), (vector * 3).components())
    print("angle: ", vector.angle(vector))

    vector = Vector.initializeFromComponents(numpy.array([11, -2, 23]))


    x1 = numpy.linspace(1, 2, 11)
    y1 = numpy.linspace(2, 3, 11)
    z1 = numpy.linspace(3, 4, 11)
    # vector = Vector(x1, y1, z1)
    vector = Vector(x1*0, x1*0, z1*0+1)
    # print("shape", vector.components().shape)
    # print("components: ", vector.components())
    #
    # vector = Vector.initializeFromComponents( [x1, y1, z1] )
    # print("shape", vector.components().shape)
    # print("components: ", vector.components())
    # print("3 v : ", vector.scalarMultiplication(3).components(), (vector * 3).components())
    #
    # print("components[0]: ", vector.components()[0])
    # print("dot: ", (vector.scalarProduct(vector)))
    #
    # print("v+v: ", (vector.addVector(vector)).components())
    # print("v+v: ", ( vector + vector).components())
    #
    # print("v/|v|: ", vector.getNormalizedVector().components())
    # print("v == v?: ", vector == vector * 2)
    # print("v == 2v?: ", vector == vector * 2)
    # print("v =", vector.toString(), vector.getX())
    #
    # print("v points=", vector.nStack(), vector.isArray())
    #
    # print("perp v", vector.getOnePerpendicularVector().components())
    #
    # print("v x v =", vector.crossProduct(vector).components())
    # print("|v| =", vector.norm())

    print("Rot(v)_{x,10deg} =", vector.rotateAroundAxis(Vector(1,0,0), -numpy.radians(10)).components())


    print("angle: ", vector.angle(vector))

