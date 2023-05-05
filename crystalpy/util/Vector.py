"""
Represents a 3d vector. Accept arrays of vectors.

vector.components() gives a numpy array with shape (3) or (3, npoints)
"""
import numpy


class Vector(object):
    def __init__(self, x, y, z):
        """
        Constructor.
        :param x: x component. It can be an array.
        :param y: y component. It can be an array of the same size of x.
        :param z: z component. It can be an array of the same size of x.
        """

        self.setComponents(x, y, z)

    @staticmethod
    def initializeFromComponents(components):
        """
        Creates a vector from a list/array of at least three elements.
        :param components: [x,y,z] components of the vector.
        :return: Vector having these x,y,z components.
        """
        return Vector(components[0],
                      components[1],
                      components[2])

    def duplicate(self):
        return Vector(self.components()[0],self.components()[1],self.components()[2])

    def setComponents(self, x, y, z):
        """
        Sets vector components.
        :param x: x component.
        :param y: y component.
        :param z: z component.
        """
        self._components = numpy.asarray([x, y, z])

    def components(self):
        """
        Returns the components of this vector a 1d three element array.
        :return:
        """
        return self._components

    def componentsStack(self):
        """
        Returns the components stack of shape (3, npoints) of this vector.
        :return:
        """
        if self.isArray():
            return self.components()
        else:
            return self.components().reshape((-1,1))


    def nStack(self):
        s = numpy.array(self.components().shape)
        if s.size == 1:
            return 1
        else:
            return s[1]

    def extractStackItem(self, i):
        x = self.getX()
        y = self.getY()
        z = self.getZ()
        return Vector.initializeFromComponents([x[i], y[i], z[i]])

    def isArray(self):
        if self.nStack() == 1:
            return False
        else:
            return True

    def getX(self):
        return self.components()[0]

    def getY(self):
        return self.components()[1]

    def getZ(self):
        return self.components()[2]

    def addVector(self, summand):
        """
        Adds two vectors.
        :param summand: The vector to add to this instance.
        :return: The sum as a vector.
        """

        wX = self.getX() + summand.getX()
        wY = self.getY() + summand.getY()
        wZ = self.getZ() + summand.getZ()
        return Vector.initializeFromComponents([wX, wY, wZ])


    def scalarMultiplication(self, k):
        """
        Scalar multiplies this vector.
        :param k: The scalar to multiply with.
        :return: Scalar multiplied vector.
        """
        return Vector.initializeFromComponents([self.getX() * k, self.getY() * k, self.getZ() * k])


    def subtractVector(self, tosubstract):
        """
        Subtract a vector from this instance.
        :param subtrahend: Vector to subtract.
        :return: The difference of the two vectors.
        """
        result = self.addVector(tosubstract.scalarMultiplication(-1.0))
        return result

    def scalarProduct(self, factor):
        """
        Calculates the scalar product of this vector with the given vector.
        :param factor: The vector to calculate the scalar product with.
        :return: Scalar product of the two vectors.
        """
        # scalar_product = numpy.dot(self.components(), factor.components())
        # scalar_product = numpy.sum( self.components() * factor.components(), axis=0)
        # return scalar_product

        wX = self.getX() * factor.getX()
        wY = self.getY() * factor.getY()
        wZ = self.getZ() * factor.getZ()
        return wX + wY + wZ


    def crossProduct(self, factor):
        """
        Calculates the cross product of two vectors.
        :param factor: The vector to form the cross product with.
        :return: Cross product of the two vectors.
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
        """
        Returns the standard norm of this norm.
        :return: Norm of this vector,
        """
        norm = self.scalarProduct(self) ** 0.5
        return norm

    def getNormalizedVector(self):
        """
        Returns a normalized vector of this vector.
        :return: Normalized vector of this vector.
        """
        return self.scalarMultiplication(self.norm() ** -1.0)

    def rotateAroundAxis(self, rotation_axis, angle):
        """
        Rotates the vector around an axis.
        :param rotation_axis: Vector specifying the rotation axis (not necessarily unit vector).
        :param angle: Rotation angle.
        :return: Rotated vector.
        """
        # For the mathematics look for: Rodrigues rotation formula.
        # http://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula
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
        """
        Returns the parallel projection of this vector along the given vector.
        :param vector: Vector defining the parallel direction.
        :return: Parallel projection along the vector.
        """
        unit_direction = vector.getNormalizedVector()
        projection_in_direction = self.scalarProduct(unit_direction)
        parallel_projection = unit_direction.scalarMultiplication(projection_in_direction)

        return parallel_projection

    def perpendicularTo(self, vector):
        """
        Returns the projection perpendicular to the given vector.
        :param vector: Vector that defines the direction.
        :return: Projection perpendicular to the given vector.
        """
        perpendicular = self.subtractVector(self.parallelTo(vector))
        return perpendicular

    def getOnePerpendicularVector(self):
        """
        Returns one arbitrary vector perpendicular to this vector.
        :return: One arbitrary vector perpendicular to this vector.
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
        """
        Return the angle between this vector and the given vector.
        :param factor: Vector to determine the angle with.
        :return: Angle between this vector and the given vector.
        """
        n1 = self.getNormalizedVector()
        n2 = factor.getNormalizedVector()

        # Determine angle between the two vectors.
        cos_angle = n1.scalarProduct(n2)
        cos_angle= numpy.clip(cos_angle, -1, 1) # just in case...
        angle = numpy.arccos(cos_angle)

        return angle

    def getVectorWithAngle(self, angle):
        """
        Returns one arbitrary vector with the given angle.
        :param angle: The requested angle.
        :return:Vector with given angle to this vector.
        """
        vector_perpendicular = self.getOnePerpendicularVector()
        vector_with_angle = self.rotateAroundAxis(vector_perpendicular, angle)

        return vector_with_angle

    def printComponents(self):
        print(self.toString())

    def toString(self):
        """
        :return: a string object containing the four components of the Stokes vector.
        """
        return "{Vx} {Vy} {Vz}".format(Vx=self.components()[0],
                                       Vy=self.components()[1],
                                       Vz=self.components()[2])


    def __eq__(self, candidate):
        """
        Determines if two vectors are equal.
        :param candidate: Vector to compare to.
        :return: True if both vectors are equal. Otherwise False.
        """
        return numpy.linalg.norm(self.components()
                              -
                              candidate.components()) < 1.e-7

    def __ne__(self, candidate):
        """
        Determines if two vectors are not equal.
        :param candidate: Vector to compare to.
        :return: True if both vectors are not equal. Otherwise False.
        """
        return not (self == candidate)

    def __add__(self, o):
        return self.addVector(o)
        # return Vector.initializeFromComponents(self.components() + o.components())

    def __sub__(self, o):
        # return Vector.initializeFromComponents(self.components() - o.components())
        return self.subtractVector(o)

    def __mul__(self, o):
        # return Vector.initializeFromComponents(self.components() * o)
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
    vector = Vector(x1, y1, z1)
    print("shape", vector.components().shape)
    print("components: ", vector.components())

    vector = Vector.initializeFromComponents( [x1, y1, z1] )
    print("shape", vector.components().shape)
    print("components: ", vector.components())
    print("3 v : ", vector.scalarMultiplication(3).components(), (vector * 3).components())

    print("components[0]: ", vector.components()[0])
    print("dot: ", (vector.scalarProduct(vector)))

    print("v+v: ", (vector.addVector(vector)).components())
    print("v+v: ", ( vector + vector).components())

    print("v/|v|: ", vector.getNormalizedVector().components())
    print("v == v?: ", vector == vector * 2)
    print("v == 2v?: ", vector == vector * 2)
    print("v =", vector.toString(), vector.getX())

    print("v points=", vector.nStack(), vector.isArray())

    print("perp v", vector.getOnePerpendicularVector().components())

    print("v x v =", vector.crossProduct(vector).components())
    print("|v| =", vector.norm())


    print("angle: ", vector.angle(vector))
