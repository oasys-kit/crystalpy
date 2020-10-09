"""
Represents a diffraction setup implementation using xraylib
Except for energy, all units are in SI.
"""

import xraylib
from crystalpy.diffraction.DiffractionSetupAbstract import DiffractionSetupAbstract

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


    def angleBragg(self, energy):
        """
        Returns the Bragg angle for a given energy.
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
                                                  self._debyeWaller, 1.0)
        return F_0

    def FH(self, energy):
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
                                                  self._debyeWaller, 1.0)
        return F_H

    def FH_bar(self, energy):
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
                                                      self._debyeWaller, 1.0)

        return F_H_bar

    def dSpacing(self):
        """
        Returns the lattice spacing d.
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
        Returns the unit cell volume.

        :return: Unit cell volume
        """
        # Retrieve unit cell volume from xraylib.
        unit_cell_volume = self._crystal['volume']

        return unit_cell_volume

if __name__ == "__main__":
    from crystalpy.diffraction.GeometryType import BraggDiffraction
    a = DiffractionSetup(geometry_type=BraggDiffraction, crystal_name="Si", thickness=1e-5,
                 miller_h=1, miller_k=1, miller_l=1,
                 asymmetry_angle=0.0,
                 azimuthal_angle=0.0,)

    print(a.angleBragg(8000.0))