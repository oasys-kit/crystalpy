"""
Represents a diffraction setup implementation using xraylib
photon energy in eV
dSpacing returns A
units are in SI.

"""
import xraylib
from crystalpy.diffraction.DiffractionSetupAbstract import DiffractionSetupAbstract
# from crystalpy.util.vector import Vector

class DiffractionSetup(DiffractionSetupAbstract):

    def __init__(self,
                 geometry_type=None, crystal_name="", thickness=1e-6,
                 miller_h=1, miller_k=1, miller_l=1,
                 asymmetry_angle=0.0,
                 azimuthal_angle=0.0,):
        """
        Constructor.
        :param geometry_type: GeometryType (BraggDiffraction,...).
        :param crystal_name: The name of the crystal, e.g. Si.
        :param thickness: The crystal thickness.
        :param miller_h: Miller index H.
        :param miller_k: Miller index K.
        :param miller_l: Miller index L.
        :param asymmetry_angle: The asymmetry angle between surface normal and Bragg normal (radians).
        :param azimuthal_angle: The angle between the projection of the Bragg normal
                                on the crystal surface plane and the x axis (radians).
        """
        super().__init__(geometry_type=geometry_type,
                         crystal_name=crystal_name,
                         thickness=thickness,
                         miller_h=miller_h, miller_k=miller_k, miller_l=miller_l,
                         asymmetry_angle=asymmetry_angle,azimuthal_angle=azimuthal_angle)

        # Load crystal from xraylib.
        self._crystal = xraylib.Crystal_GetCrystal(self.crystalName())


    # todo: allow energy be a vector (like in DiffractionSetupDabax) ???
    def angleBragg(self, energy):
        """
        Returns the Bragg angle in rad for a given energy.
        :param energy: Energy to calculate the Bragg angle for.
        :return: Bragg angle.
        """
        energy_in_kev = energy / 1000.0

        # Retrieve bragg angle from xraylib.
        angle_bragg = xraylib.Bragg_angle(self._crystal,
                                          energy_in_kev,
                                          self.millerH(),
                                          self.millerK(),
                                          self.millerL())
        return angle_bragg

    def F0(self, energy):
        """
        Calculate F0 from Zachariasen.
        :param energy: photon energy in eV.
        :return: F0
        """
        energy_in_kev = energy / 1000.0
        F_0 = xraylib.Crystal_F_H_StructureFactor(self._crystal,
                                                  energy_in_kev,
                                                  0, 0, 0,
                                                  self._debyeWaller, rel_angle=0.0)
        return F_0

    def FH(self, energy, rel_angle=1.0):
        """
        Calculate FH from Zachariasen.
        :param energy: photon energy in eV.
        :return: FH
        """
        energy_in_kev = energy / 1000.0
        F_H = xraylib.Crystal_F_H_StructureFactor(self._crystal,
                                                  energy_in_kev,
                                                  self.millerH(),
                                                  self.millerK(),
                                                  self.millerL(),
                                                  self._debyeWaller, rel_angle)
        return F_H

    def FH_bar(self, energy, rel_angle=1.0):
        """
        Calculate FH_bar from Zachariasen.
        :param energy: photon energy in eV.
        :return: FH_bar
        """
        energy_in_kev = energy / 1000.0
        F_H_bar = xraylib.Crystal_F_H_StructureFactor(self._crystal,
                                                      energy_in_kev,
                                                      -self.millerH(),
                                                      -self.millerK(),
                                                      -self.millerL(),
                                                      self._debyeWaller, rel_angle)

        return F_H_bar

    def Fall(self, energy, rel_angle=1.0):
        return self.F0(energy),\
               self.FH(energy, rel_angle=rel_angle), \
               self.FH_bar(energy, rel_angle=rel_angle)

    def dSpacing(self):
        """
        Returns the lattice spacing d in A.
        :return: Lattice spacing.
        """

        # Retrieve lattice spacing d from xraylib in Angstrom.
        d_spacing = xraylib.Crystal_dSpacing(self._crystal,
                                             self.millerH(),
                                             self.millerK(),
                                             self.millerL())

        return d_spacing

    def unitcellVolume(self):
        """
        Returns the unit cell volume in A**3.

        :return: Unit cell volume
        """
        # Retrieve unit cell volume from xraylib.
        unit_cell_volume = self._crystal['volume']

        return unit_cell_volume


if __name__ == "__main__":
    from crystalpy.diffraction.GeometryType import BraggDiffraction
    import numpy

    a = DiffractionSetup(geometry_type=BraggDiffraction, crystal_name="Si", thickness=1e-5,
                 miller_h=1, miller_k=1, miller_l=1,
                 asymmetry_angle=0.0,
                 azimuthal_angle=0.0,)

    energy = 8000.0
    print("Photon energy: %g deg " % (energy))
    print("d_spacing: %g A " % (a.dSpacing()))
    print("unitCellVolumw: %g A**3 " % (a.unitcellVolume()))
    print("F0 ", a.F0(energy))
    print("FH ", a.FH(energy))
    print("FH_BAR ", a.FH_bar(energy))

    print("PSI0 ", a.psi0(energy))
    print("PSIH ", a.psiH(energy))
    print("PSIH_bar ", a.psiH_bar(energy))

    print("V0: ", a.vectorK0direction(energy).components())
    print("Bh direction: ", a.vectorHdirection().components())
    print("Bh: ", a.vectorH().components())
    print("K0: ", a.vectorK0(energy).components())
    print("Kh: ", a.vectorKh(energy).components())
    print("Vh: ", a.vectorKhdirection(energy).components())


    from crystalpy.util.Photon import Photon
    print("Difference to ThetaB uncorrected: ",
          a.deviationOfIncomingPhoton(Photon(energy_in_ev=energy, direction_vector=a.vectorK0(energy))))


    # print("Asymmerey factor b: ", a.asymmetry_factor(energy))
    print("Bragg angle: %g deg " %  (a.angleBragg(energy) * 180 / numpy.pi))
    print("Bragg angle corrected: %g deg " %  (a.angleBraggCorrected(energy) * 180 / numpy.pi))


 #     VIN_BRAGG_UNCORR (Uncorrected): (  0.00000000,    0.968979,   -0.247145)
 #     VIN_BRAGG          (Corrected): (  0.00000000,    0.968971,   -0.247176)
 #     VIN_BRAGG_ENERGY              : (  0.00000000,    0.968971,   -0.247176)
 # Reflected directions matching Bragg angle:
 #    VOUT_BRAGG_UNCORR (Uncorrected): (  0.00000000,    0.968979,    0.247145)
 #    VOUT_BRAGG          (Corrected): (  0.00000000,    0.968971,    0.247176)
 #    VOUT_BRAGG_ENERGY              : (  0.00000000,    0.968971,    0.247176)
