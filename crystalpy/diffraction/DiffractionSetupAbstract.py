"""
Represents a diffraction setup abstract class.
The super class should implement the methods to calculate structure factors
"""

from collections import OrderedDict
from copy import deepcopy
import numpy
import scipy.constants as codata

from crystalpy.util.Vector import Vector


class DiffractionSetupAbstract(object):
    """
    Constructor.

    Parameters
    ----------
    geometry_type: instance of BraggDiffraction, LaueDiffraction, BraggTransmission, or LaueTransmission

    crystal_name: str
        The name of the crystal, e.g. "Si".

    thickness: float
        The crystal thickness in m.

    miller_h: int
        Miller index H.

    miller_k: int
        Miller index K.

    miller_l: int
        Miller index L.

    asymmetry_angle: float
        The asymmetry angle between surface normal and Bragg normal (radians).

    azimuthal_angle: float
        The angle between the projection of the Bragg normal on the crystal surface plane and the x axis (radians).

    debye_waller: float
        The Debye-Waller factor exp(-M).

    """
    def __init__(self,
                 geometry_type=None,
                 crystal_name="",
                 thickness=1e-6,
                 miller_h=1,
                 miller_k=1,
                 miller_l=1,
                 asymmetry_angle=0.0,
                 azimuthal_angle=0.0,
                 debye_waller=1.0):

        self._geometry_type = geometry_type
        self._crystal_name = crystal_name
        self._thickness = thickness
        self._miller_h = miller_h
        self._miller_k = miller_k
        self._miller_l = miller_l
        self._asymmetry_angle = asymmetry_angle  # degrees.
        self._azimuthal_angle = azimuthal_angle  # degrees
        self._debyeWaller = debye_waller

    #
    # setters and getters
    #
    def geometryType(self):
        """Returns the GeometryType, e.g. BraggDiffraction, LaueTransmission,...

        Returns
        -------
        instance of BraggDiffraction, LaueDiffraction, BraggTransmission, or LaueTransmission.
            The GeometryType.


        """
        return self._geometry_type

    def crystalName(self):
        """Returs the crystal name

        Returns
        -------
        str
            Crystal name.

        """
        return self._crystal_name

    def thickness(self):
        """Returns the crystal thickness in meters

        Returns
        -------
        float
            he crystal thickness.

        """
        return self._thickness

    def millerH(self):
        """Returns the Miller H index.

        Returns
        -------
        int
            Miller H index.

        """
        return self._miller_h

    def millerK(self):
        """Returns the Miller K index.

        Returns
        -------
        int
            Miller K index.

        """
        return self._miller_k

    def millerL(self):
        """Returns the Miller L index.

        Returns
        -------
        int
            Miller L index.

        """
        return self._miller_l

    def asymmetryAngle(self):
        """Returns the asymmetry angle between surface normal and Bragg normal in degrees.

        Returns
        -------
        float
            Asymmetry angle.


        """
        return self._asymmetry_angle

    def azimuthalAngle(self):
        """Returns the angle between the Bragg normal projection on the crystal surface plane and the x axis in degrees.

        Returns
        -------
        float
            Azimuthal angle.

        """
        return self._azimuthal_angle

    #
    # abstract methods to be implemented by the super class
    #
    def angleBragg(self, energy=8000.0):
        """Returns the Bragg angle for a given energy in radians.

        Parameters
        ----------
        energy :
            Energy to calculate the Bragg angle for. (Default value = 8000.0)

        Returns
        -------
        float
            Bragg angle in radians.

        """
        raise NotImplementedError()

    def F0(self, energy=8000.0):
        """Calculate the structure factor F0.

        Parameters
        ----------
        energy : float
            photon energy in eV. (Default value = 8000.0)

        Returns
        -------
        complex
            F0

        """
        raise NotImplementedError()


    def FH(self, energy=8000.0):
        """Calculate the structure factor FH.

        Parameters
        ----------
        energy :
            photon energy in eV. (Default value = 8000.0)

        Returns
        -------
        complex
            FH

        """

        raise NotImplementedError()

    def FH_bar(self, energy=8000.0):
        """Calculate the structure factor  FH_bar.

        Parameters
        ----------
        energy :
            photon energy in eV. (Default value = 8000.0)

        Returns
        -------
        complex
            FH_bar

        """

        raise NotImplementedError()

    def Fall(self, energy=8000.0):
        """Calculate the all structure factor  (F0, FH, FH_bar).

        Parameters
        ----------
        energy :
            photon energy in eV. (Default value = 8000.0)

        Returns
        -------
        tuple
            (F0, FH, FH_bar).

        """

        raise NotImplementedError()

    def dSpacing(self):
        """Returns the lattice spacing d in A.

        Returns
        -------
        float
            Lattice spacing. in A

        """

        raise NotImplementedError()


    def unitcellVolume(self):
        """Returns the unit cell volume in A^3

        Returns
        -------
        float
            Unit cell volume in A^3.

        """
        raise NotImplementedError()

    #
    # other methods
    #

    def dSpacingSI(self):
        """Returns the lattice spacing d in SI units (meters).

        Returns
        -------
        float
            Lattice spacing. in m

        """
        return 1e-10 * self.dSpacing()

    def unitcellVolumeSI(self):
        """Returns the unit cell volume in SI units (m^3)

        Returns
        -------
        float
            Unit cell volume in m^3.

        """
        return 1e-30 * self.unitcellVolume()
    #
    # structure factors
    #
    def psi0(self, energy):
        """Calculate the structure factor psi0 (defined in Zachariasen [3-95]).

        Parameters
        ----------
        energy :
            photon energy in eV. (Default value = 8000.0)

        Returns
        -------
        complex
            psi0

        """
        classical_electron_radius = codata.codata.physical_constants["classical electron radius"][0]
        wavelength = codata.h * codata.c / codata.e / energy
        return (-classical_electron_radius * wavelength ** 2 / (numpy.pi * self.unitcellVolumeSI())) * self.F0(energy)

    def psiH(self, energy, rel_angle=1.0):
        """Calculate the structure factor psiH (defined in Zachariasen [3-95]).

        Parameters
        ----------
        energy :
            photon energy in eV. (Default value = 8000.0)

        rel_angle : float, optional
            (Default = 1.0)

        Returns
        -------
        complex
            psiH

        """
        classical_electron_radius = codata.codata.physical_constants["classical electron radius"][0]
        wavelength = codata.h * codata.c / codata.e / energy
        return (-classical_electron_radius * wavelength ** 2 / (numpy.pi * self.unitcellVolumeSI())) * self.FH(energy, rel_angle=rel_angle)

    def psiH_bar(self, energy, rel_angle=1.0):
        """Calculate the structure factor psiH_bar (defined in Zachariasen [3-95]).

        Parameters
        ----------
        energy :
            photon energy in eV. (Default value = 8000.0)

        rel_angle : float, optional
            (Default = 1.0)

        Returns
        -------
        complex
            psiH_bar

        """
        classical_electron_radius = codata.codata.physical_constants["classical electron radius"][0]
        wavelength = codata.h * codata.c / codata.e / energy
        return (-classical_electron_radius * wavelength ** 2 / (numpy.pi * self.unitcellVolumeSI())) * self.FH_bar(energy, rel_angle=rel_angle)

    def psiAll(self, energy1, rel_angle=1.0):
        """Calculate the psi structure factors (psi0, psiH, psiH_bar) (defined in Zachariasen [3-95]).

        Parameters
        ----------
        energy :
            photon energy in eV. (Default value = 8000.0)

        rel_angle : float, optional
            (Default = 1.0)

        Returns
        -------
        tuple
            (psi0, psiH, psiH_bar).

        """
        energy = numpy.array(energy1)
        classical_electron_radius = codata.codata.physical_constants["classical electron radius"][0]
        wavelength = codata.h * codata.c / codata.e / energy
        factor = (-classical_electron_radius * wavelength ** 2 / (numpy.pi * self.unitcellVolumeSI()))
        Fall = self.Fall(energy, rel_angle=rel_angle)
        return  factor*Fall[0], factor*Fall[1], factor*Fall[2]

    #
    # vector interface
    #
    def vectorNormalSurface(self):
        """Returns the normal to the surface. (0,0,1) by definition.

        Returns
        -------
        Vector instance
            Vector instance with Surface normal Vnor.

        """
        # Geometrical convention from M.Sanchez del Rio et al., J.Appl.Cryst.(2015). 48, 477-491.
        return Vector(0, 0, 1)

    def vectorNormalSurfaceInwards(self):
        """Returns the inwards normal to the surface. -vectorNormalSurface() by definition.

        Returns
        -------
        Vector instance
            Vector instance with inwards surface normal Vnor.

        """
        # Geometrical convention from M.Sanchez del Rio et al., J.Appl.Cryst.(2015). 48, 477-491.
        return self.vectorNormalSurface().scalarMultiplication(-1.0)

    def vectorParallelSurface(self):
        """Returns the direction parallel to the crystal surface. (0,1,0) by definition.

        Returns
        -------
        vector instance
            Vector instance with Surface normal Vtan.

        """
        # Geometrical convention from M.Sanchez del Rio et al., J.Appl.Cryst.(2015). 48, 477-491.
        return  Vector(0, 1, 0)

    def vectorH(self):
        """Calculates the H vector, normal on the reflection lattice plane, with modulus 2 pi / d_spacing (SI).
        
        The normal to Bragg planes is obtained by rotating vnor an angle equal to minuns asymmetry angle (-alphaXOP)
        around X using rodrigues rotation (in the screw direction (cw) when looking in the axis direction),
        and then an angle phi (azimuthal angle) around Z


        Returns
        -------
        vector instance
            H vector

        References
        ----------
        Sanchez del Rio, M., Perez-Bocanegra, N., Shi, X., Honkimäki, V. & Zhang, L. (2015).
        Simulation of X-ray diffraction profiles for bent anisotropic crystals. J. Appl. Cryst. 48, 477–491.
        http://dx.doi.org/10.1107/S1600576715002782

        """
        # Geometrical convention from M.Sanchez del Rio et al., J.Appl.Cryst.(2015). 48, 477-491.

        g_modulus = 2.0 * numpy.pi / (self.dSpacingSI())
        # Let's start from a vector parallel to the surface normal (z axis).
        temp_normal_bragg = Vector(0, 0, 1).scalarMultiplication(g_modulus)

        # Let's now rotate this vector of an angle alphaX around the y axis (according to the right-hand-rule).
        alpha_x = self.asymmetryAngle()
        axis = self.vectorParallelSurface().crossProduct(self.vectorNormalSurface())  # should be Vector(1, 0, 0)
        temp_normal_bragg = temp_normal_bragg.rotateAroundAxis(axis, -alpha_x)

        # Let's now rotate this vector of an angle phi around the z axis (following the ISO standard 80000-2:2009).
        phi = self.azimuthalAngle()
        normal_bragg = temp_normal_bragg.rotateAroundAxis(Vector(0, 0, 1), phi)

        return normal_bragg

    def vectorHdirection(self):
        """Calculates the unitary vector parallel to the H vector (normal on the reflection lattice plane, with modulus 2 pi / d_spacing (SI)).

        The normal to the Bragg planes is obtained by rotating vnor an angle equal to minuns asymmetry angle (-alphaXOP)
        around X using rodrigues rotation (in the screw direction (cw) when looking in the axis direction),
        and then an angle phi (azimuthal angle) around Z

        Returns
        -------
        Vector instance
            normal vector in direction of H.

        """
        return self.vectorH().getNormalizedVector()

    def vectorK0direction(self, energy):
        """Calculates the unitary vector parallel to the K0 vector (along the Bragg position)

        Parameters
        ----------
        energy : float or numpy array.
            The photon energy in eV.
            

        Returns
        -------
        Vector instance
            The normalized vector (or stack of vectors) with the directions of K0.

        """
        minusBH = self.vectorHdirection().scalarMultiplication(-1.0) # -BH of an angle (90-BraggAngle) around the x axis
        axis = self.vectorParallelSurface().crossProduct(self.vectorNormalSurface())  # should be Vector(1, 0, 0)
        photon_direction = minusBH.rotateAroundAxis(axis, (numpy.pi / 2) - self.angleBragg(energy))
        return photon_direction

    def vectorK0directionCorrected(self, energy):
        """Calculates the unitary vector parallel to the K0corrected vector (along the Bragg position corrected for refraction)

        Parameters
        ----------
        energy : float or numpy array.
            The photon energy in eV.


        Returns
        -------
        Vector instance
            The normalized vector (or stack of vectors) with the directions of K0corrected.

        """
        minusBH = self.vectorHdirection().scalarMultiplication(-1.0) # -BH of an angle (90-BraggAngle) around the x axis
        axis = self.vectorParallelSurface().crossProduct(self.vectorNormalSurface())  # should be Vector(1, 0, 0)
        photon_direction = minusBH.rotateAroundAxis(axis, (numpy.pi / 2) - self.angleBraggCorrected(energy))
        return photon_direction

    def vectorK0(self, energy):
        """Calculates the vector K0(along the Bragg position)

        Parameters
        ----------
        energy : float or numpy array.
            The photon energy in eV.


        Returns
        -------
        Vector instance
            The K0.

        """
        wavelength = codata.h * codata.c / codata.e / energy
        return self.vectorK0direction(energy).scalarMultiplication(2 * numpy.pi / wavelength)

    def vectorK0corrected(self, energy):
        """Calculates the vector K0corrected (along the corrected Bragg position)

        Parameters
        ----------
        energy : float or numpy array.
            The photon energy in eV.


        Returns
        -------
        Vector instance
            The K0corrected.

        """
        wavelength = codata.h * codata.c / codata.e / energy
        return self.vectorK0directionCorrected(energy).scalarMultiplication(2 * numpy.pi / wavelength)

    def vectorKh(self, energy):
        """returns KH that verifies Laue equation with K0

        Parameters
        ----------
        energy : float or numpy array
            The energy or energy array

        Returns
        -------
        Vector instance
            The KH vector or vector stack

        """
        return Vector.addVector(self.vectorK0(energy), self.vectorH())

    def vectorKhdirection(self, energy):
        """returns an unitary vector along the KH direction (that that verifies Laue equation with K0).

        Parameters
        ----------
        energy : float or numpy array
            The energy or energy array

        Returns
        -------
        Vector instance
            The unitary vector(s) along the KH direction(s).

        """
        return self.vectorKh(energy).getNormalizedVector()

    def vectorKscattered(self, K_IN=None, energy=8000.0):
        """
        returns the scattered K vector following the scattering equation at a surface:
            K_parallel = K_IN_parallel + H_parallel
            |K| = |K_IN|

        Parameters
        ----------
        K_IN : instance of Vector, optional
            The K vector. If None, used the vectorK0corrected(energy)
        energy : float, optional
            The energy value in eV (used only if K_IN=None)

        Returns
        -------
        Vector instance
            Vector with the scattered K.

        """
        if K_IN is None:
            K_IN = self.vectorK0corrected(energy)

        H = self.vectorH()
        NORMAL = self.vectorNormalSurface()
        K_OUT = K_IN.scatteringOnSurface(NORMAL, H)
        return K_OUT

    # useful for scans...
    def vectorIncomingPhotonDirection(self, energy, deviation, angle_center_flag=2):
        """Calculates the direction of the incoming photon (or photon stack). Parallel to k_0.

        Parameters
        ----------
        energy : float of numpy array/
            Energy in eV.

        deviation : float or array.
            Deviation from the uncorrected Bragg angle.
            A positive deviation means the photon direction lies closer to the surface normal.

        angle_center_flag : int, optional
             Flag from where "deviation: is measured:
             0: absolute angle, 1: from Bragg angle corrected for refraction, 2: from Bragg angle.

        Returns
        -------
        Vector instance
            Direction(s) of the incoming photon(s).

        """
        # Geometrical convention from M.Sanchez del Rio et al., J.Appl.Cryst.(2015). 48, 477-491.

        # # DONE: vectorize this part as in https://github.com/srio/CRYSTAL/blob/master/crystal3.F90
        # angle between the incoming photon direction and the surface normal (z axis).
        # a positive deviation means the photon direction lies closer to the surface normal.
        # angle = numpy.pi / 2.0 - (self.angleBragg(energy) + self.asymmetryAngle() + deviation)
        # # the photon comes from left to right in the yz plane.
        # photon_direction_old = Vector(0,numpy.sin(angle),-numpy.cos(angle))
        # angle_center_flag = 0,  # 0=Absolute angle, 1=Theta Bragg Corrected, 2=Theta Bragg

        # print(">>>>> in vectorIncomingPhotonDirection")

        # Let's now rotate -BH of an angle (90-BraggAngle) around the x axis
        minusBH = self.vectorH().scalarMultiplication(-1.0)
        minusBH = minusBH.getNormalizedVector()
        axis = self.vectorParallelSurface().crossProduct(self.vectorNormalSurface())  # should be Vector(1, 0, 0)

        if angle_center_flag == 0:
            photon_direction = minusBH.rotateAroundAxis(axis, (numpy.pi / 2) - deviation)
        elif angle_center_flag == 1:
            photon_direction = minusBH.rotateAroundAxis(axis, (numpy.pi / 2) - self.angleBraggCorrected(energy) - deviation)
        elif angle_center_flag == 2:
            photon_direction = minusBH.rotateAroundAxis(axis, (numpy.pi / 2) - self.angleBragg(energy) - deviation)

        # print("PHOTON DIRECTION ",photon_direction_old.components(),photon_direction.components())
        # Let's now rotate this vector of an angle phi around the z axis (following the ISO standard 80000-2:2009).
        # photon_direction = photon_direction.rotateAroundAxis(Vector(0, 0, 1), self.azimuthalAngle() )

        return photon_direction

    # def vectorIncomingPhotonDirection(self, energy, deviation):
    #     """Calculates the direction of the incoming photon. Parallel to K0.
    #
    #     Parameters
    #     ----------
    #     energy : float or numpy array
    #         Energy in eV.
    #
    #     deviation : float or numpy array
    #         Deviation from the Bragg angle in radians.
    #
    #     Returns
    #     -------
    #     Vector instance
    #         Direction(s) of the incoming photon(s).
    #
    #     """
    #     # Edoardo: I use the geometrical convention from
    #     # M.Sanchez del Rio et al., J.Appl.Cryst.(2015). 48, 477-491.
    #
    #     # # DONE: vectorize this part as in https://github.com/srio/CRYSTAL/blob/master/crystal3.F90
    #     # # angle between the incoming photon direction and the surface normal (z axis).
    #     # # a positive deviation means the photon direction lies closer to the surface normal.
    #     # angle = numpy.pi / 2.0 - (self.angleBragg(energy) + self.asymmetryAngle() + deviation)
    #     # # the photon comes from left to right in the yz plane.
    #     # photon_direction_old = Vector(0,numpy.sin(angle),-numpy.cos(angle))
    #
    #     print(">>>>> ****** in vectorIncomingPhotonDirection")
    #     # Let's now rotate -BH of an angle (90-BraggAngle) around the x axis
    #     minusBH = self.vectorHdirection().scalarMultiplication(-1.0)
    #     # minusBH = minusBH.getNormalizedVector()
    #     axis = self.vectorParallelSurface().crossProduct(self.vectorNormalSurface())  # should be Vector(1, 0, 0)
    #     # TODO check why deviation has minus
    #     photon_direction = minusBH.rotateAroundAxis(axis, (numpy.pi/2)-self.angleBragg(energy)-deviation)
    #
    #     # print("PHOTON DIRECTION ",photon_direction_old.components(),photon_direction.components())
    #     # Let's now rotate this vector of an angle phi around the z axis (following the ISO standard 80000-2:2009).
    #     # photon_direction = photon_direction.rotateAroundAxis(Vector(0, 0, 1), self.azimuthalAngle() )
    #
    #     return photon_direction

    #
    # tools
    #
    def clone(self):
        """Returns a copy of this instance.

        Returns
        -------
        DiffractionSetup instance
            A copy of this instance.

        """
        return deepcopy(self)

    def duplicate(self):
        """Returns a copy of this instance.

        Returns
        -------
        DiffractionSetup instance
            A copy of this instance.

        """
        return deepcopy(self)

    def toDictionary(self):
        """Returns info of this setup in a dictionary.

        Returns
        -------
        dict
            Info dictionary form of this setup.

        """
        info_dict = OrderedDict()
        info_dict["Geometry Type"] = self.geometryType().description()
        info_dict["Crystal Name"] = self.crystalName()
        info_dict["Thickness"] = str(self.thickness())
        info_dict["Miller indices (h,k,l)"] = "(%i,%i,%i)" % (self.millerH(),
                                                              self.millerK(),
                                                              self.millerL())
        info_dict["Asymmetry Angle"] = str(self.asymmetryAngle())
        info_dict["Azimuthal Angle"] = str(self.azimuthalAngle())

        return info_dict


    def deviationOfIncomingPhoton(self, photon_in):
        """Calculates deviation from the Bragg angle of an incoming photon in radians.

        Parameters
        ----------
        photon_in :
            Incoming photon.

        Returns
        -------
        float
            Deviation from Bragg angle in radians.

        """
        # this holds for every incoming photon-surface normal plane.
        total_angle = photon_in.unitDirectionVector().angle(self.vectorH())

        energy = photon_in.energy()
        angle_bragg = self.angleBragg(energy)

        deviation = total_angle - angle_bragg - numpy.pi / 2
        return deviation


    # """
    # ! asymmetry b factor vectorial value (Zachariasen, [3.115])
    # """

    def asymmetryFactor(self, energy, vector_k_in=None):
        """Returns asymmetric factor (after Zachariasen equation [3.115]).

        Parameters
        ----------
        energy : float or numpy array
            The photon energy in eV.
            
        vector_k_in : Vector instance, optional
             The incident K0 (Default value = None, meaning that K0 is used.)

        Returns
        -------
        float or numpy array

        """
        if vector_k_in is None:
            vector_k_in = self.vectorK0(energy)

        v2 = vector_k_in.addVector(self.vectorKh(energy)).subtractVector(self.vectorK0(energy))

        numerator = Vector.scalarProduct(self.vectorNormalSurfaceInwards(),vector_k_in)
        denominator = Vector.scalarProduct(self.vectorNormalSurfaceInwards(),v2)

        return numerator / denominator

    def angleBraggCorrected(self, energy=8000.0, use_exact_equation=True):
        """Returns the Bragg angle corrected for refraction for a given energy.
        An approximated formula is found in Zachariasen equation 3.145a.
        The exact formula is in Guigay % Sanchez del Rio equation 21.

        Parameters
        ----------
        energy : float or numpy array
            Energy in eV for calculating the Bragg angle. (Default value = 8000.0)

        Returns
        -------
        float or numpy array
            Bragg angle(s) corrected.

        """

        if use_exact_equation:
            numerator = (1 - self.asymmetryFactor(energy)) * self.psi0(energy).real
            denominator = 4 * self.asymmetryFactor(energy) * numpy.sin(self.angleBragg(energy))
            # equation 21 in G&SR
            return numpy.arcsin( numpy.sin(self.angleBragg(energy)) + numerator / denominator)
        else:
            numerator = (1 - self.asymmetryFactor(energy)) * self.psi0(energy).real
            denominator = 2 * self.asymmetryFactor(energy) * numpy.sin(2 * self.angleBragg(energy))
            # equation 3.145a in Zachariasen's book
            return self.angleBragg(energy) + numerator / denominator

    #
    # Darwin width
    #

    def darwinHalfwidthS(self, energy=8000.0):
        """

        energy : float or numpy array
            Energy in eV for calculating the Bragg angle. (Default value = 8000.0)
            

        Returns
        -------
        float or numpy array
            1/2 of the Darwin width(s) for sigma polarization in radians.

        """
        return self.darwinHalfwidth(energy)[0]

    def darwinHalfwidthP(self, energy):
        """

        energy : float or numpy array
            Energy in eV for calculating the Bragg angle. (Default value = 8000.0)


        Returns
        -------
        float or numpy array
            1/2 of the Darwin width(s) for pi polarization in radians.

        """
        return self.darwinHalfwidth(energy)[1]

    def darwinHalfwidth(self, energy):
        """

        Parameters
        ----------
        energy :
            

        Returns
        -------

        """
        if isinstance(energy, int): energy = float(energy)

        codata_e2_mc2 = codata.hbar * codata.alpha / codata.m_e / codata.c * 1e2 # in cm
        wavelength = codata.c * codata.h / codata.e / energy

        RN = 1.0 / (self.unitcellVolumeSI() * 1e6 ) * codata_e2_mc2
        R_LAM0 = wavelength * 1e2
        F_0, FH, FH_BAR = self.Fall(energy)
        STRUCT = numpy.sqrt( FH * FH_BAR)
        TEMPER = 1.0 # self.get_preprocessor_dictionary()["temper"]
        GRAZE = self.angleBragg(energy)
        SSVAR	= RN * (R_LAM0**2) * STRUCT * TEMPER / numpy.pi / numpy.sin(2.0 * GRAZE)
        SPVAR = SSVAR * numpy.abs(numpy.cos(2.0 * GRAZE))
        return SSVAR.real, SPVAR.real

    #
    # operators
    #
    def __eq__(self, candidate):
        if self._geometry_type != candidate.geometryType():
            return False

        if self._crystal_name != candidate.crystalName():
            return False

        if self._thickness != candidate.thickness():
            return False

        if self._miller_h != candidate.millerH():
            return False

        if self._miller_k != candidate.millerK():
            return False

        if self._miller_l != candidate.millerL():
            return False

        if self._asymmetry_angle != candidate.asymmetryAngle():
            return False

        if self._azimuthal_angle != candidate.azimuthalAngle():
            return False

        # All members are equal so are the instances.
        return True

    def __ne__(self, candidate):
        return not self == candidate



if __name__ == "__main__":
    if False:
        a = DiffractionSetupAbstract(geometry_type=0, crystal_name="Si", thickness=1e-5,
                     miller_h=1, miller_k=1, miller_l=1,
                     asymmetry_angle=0.0,
                     azimuthal_angle=0.0,)
