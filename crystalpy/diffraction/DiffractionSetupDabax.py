"""
Represents a diffraction setup implementation using DABAX
photon energy in eV
dSpacing returns A
units are in SI.
"""

from crystalpy.diffraction.DiffractionSetupAbstract import DiffractionSetupAbstract
from dabax.dabax_xraylib import DabaxXraylib

class DiffractionSetupDabax(DiffractionSetupAbstract):
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
        The angle between the projection of the Bragg normal on the crystal surface plane and the Y axis (radians).
        It can also be called inclination angle.

    debye_waller: float
        The Debye-Waller factor exp(-M).

    crystal_data: None or dict
        If None, the crystal data is loaded from dabax Crystals.dat file (entry: crystal_name).
        Alternatively, we can force using a user-structure entered here as a dictionary with the xraylib format.

    dabax: DabaxXraylib instance, optional
        Default : None, use DabaxXraylib()

    """

    def __init__(self,
                 geometry_type=None, crystal_name="", thickness=1e-6,
                 miller_h=1, miller_k=1, miller_l=1,
                 asymmetry_angle=0.0,
                 azimuthal_angle=0.0,
                 crystal_data=None,
                 dabax=None):

        super().__init__(geometry_type=geometry_type,
                         crystal_name=crystal_name,
                         thickness=thickness,
                         miller_h=miller_h,
                         miller_k=miller_k,
                         miller_l=miller_l,
                         asymmetry_angle=asymmetry_angle,
                         azimuthal_angle=azimuthal_angle)

        if isinstance(dabax, DabaxXraylib):
            self.dx = dabax
        else:
            self.dx = DabaxXraylib()

        # Load crystal from dabax (if not defined in crystal_data).
        if crystal_data is None:
            self._crystal = self.dx.Crystal_GetCrystal(self.crystalName())
        else:
            self._crystal = crystal_data

    def angleBragg(self, energy):
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
        energy_in_kev = energy / 1000.0


        angle_bragg = self.dx.Bragg_angle(self._crystal,
                                          energy_in_kev,
                                          self.millerH(),
                                          self.millerK(),
                                          self.millerL(),)

        return angle_bragg

    def F0(self, energy):
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
        return self.Fall(energy)[0]


    def FH(self, energy, rel_angle=1.0):
        """Calculate the structure factor FH.

        Parameters
        ----------
        energy :
            photon energy in eV. (Default value = 8000.0)

        rel_angle: float, optional
            ratio of the incident angle and the Bragg angle (Default : rel_angle=1.0)

        Returns
        -------
        complex
            FH

        """
        return self.Fall(energy, rel_angle=rel_angle)[1]


    def FH_bar(self, energy, rel_angle=1.0):
        """Calculate the structure factor  FH_bar.

        Parameters
        ----------
        energy :
            photon energy in eV. (Default value = 8000.0)

        rel_angle: float, optional
            ratio of the incident angle and the Bragg angle (Default : rel_angle=1.0)

        Returns
        -------
        complex
            FH_bar

        """
        return self.Fall(energy, rel_angle=rel_angle)[2]


    def Fall(self, energy, rel_angle=1.0):
        """Calculate the all structure factor  (F0, FH, FH_bar).

        Parameters
        ----------
        energy :
            photon energy in eV. (Default value = 8000.0)

        rel_angle: float, optional
            ratio of the incident angle and the Bragg angle (Default : rel_angle=1.0)

        Returns
        -------
        tuple
            (F0, FH, FH_bar).

        """

        energy_in_kev = energy / 1000.0

        Fall = self.dx.Crystal_F_0_F_H_F_H_bar_StructureFactor(self._crystal,
                                                  energy_in_kev,
                                                  self.millerH(),
                                                  self.millerK(),
                                                  self.millerL(),
                                                  self._debyeWaller, rel_angle)
        return Fall


    def dSpacing(self):
        """Returns the lattice spacing d in A.

        Returns
        -------
        float
            Lattice spacing. in A

        """

        # Retrieve lattice spacing d from xraylib in Angstrom.

        d_spacing = self.dx.Crystal_dSpacing(self._crystal,
                                             self.millerH(),
                                             self.millerK(),
                                             self.millerL())


        return d_spacing

    def unitcellVolume(self):
        """Returns the unit cell volume in A^3

        Returns
        -------
        float
            Unit cell volume in A^3.

        """
        # Retrieve unit cell volume from xraylib.
        unit_cell_volume = self._crystal['volume']

        return unit_cell_volume


if __name__ == "__main__":

    if False:
        from crystalpy.diffraction.GeometryType import BraggDiffraction
        import numpy

        a = DiffractionSetupDabax(geometry_type=BraggDiffraction, crystal_name="Si", thickness=1e-5,
                     miller_h=1, miller_k=1, miller_l=1,
                     asymmetry_angle=0.0,
                     azimuthal_angle=0.0,
                     dabax=DabaxXraylib())

        energy = 8000.0
        energies = numpy.linspace(energy, energy + 100, 2)
        print("Photon energy: %g deg " % (energy))
        print("d_spacing: %g A " % (a.dSpacing()))
        print("unitCellVolumw: %g A**3 " % (a.unitcellVolume()))
        print("Bragg angle: %g deg " %  (numpy.degrees(a.angleBragg(energy))))
        print("Bragg angle [array] [deg] ", numpy.degrees(a.angleBragg(energies)))
        print("Asymmerey factor b: ", a.asymmetryFactor(energy))

        print("F0 ", a.F0(energy))
        print("F0 [array] ", a.F0(energies))
        print("FH ", a.FH(energy))
        print("FH [array] ", a.FH(energies))
        print("FH_bar ", a.FH_bar(energy))
        print("FH_bar [array] ", a.FH_bar(energies))

        print("PSI0 ", a.psi0(energy))
        print("PSIH ", a.psiH(energy))
        print("PSIH_bar ", a.psiH_bar(energy))
        #
        print("V0: ", a.vectorK0direction(energy).components())
        print("Bh direction: ", a.vectorHdirection().components())
        print("Bh: ", a.vectorH().components())
        print("K0: ", a.vectorK0(energy).components())
        print("Kh: ", a.vectorKh(energy).components())
        print("Vh: ", a.vectorKhdirection(energy).components())
        #
        #
        from crystalpy.util.Photon import Photon
        print("Difference to ThetaB uncorrected: ",
              a.deviationOfIncomingPhoton(Photon(energy_in_ev=energy, direction_vector=a.vectorK0(energy))))
        #
        #
        print("Asymmerey factor b: ", a.asymmetryFactor(energy))
        print("Bragg angle: %g deg " %  (a.angleBragg(energy) * 180 / numpy.pi))
        # print("Bragg angle corrected: %g deg " %  (a.angleBraggCorrected(energy) * 180 / numpy.pi))


     #     VIN_BRAGG_UNCORR (Uncorrected): (  0.00000000,    0.968979,   -0.247145)
     #     VIN_BRAGG          (Corrected): (  0.00000000,    0.968971,   -0.247176)
     #     VIN_BRAGG_ENERGY              : (  0.00000000,    0.968971,   -0.247176)
     # Reflected directions matching Bragg angle:
     #    VOUT_BRAGG_UNCORR (Uncorrected): (  0.00000000,    0.968979,    0.247145)
     #    VOUT_BRAGG          (Corrected): (  0.00000000,    0.968971,    0.247176)
     #    VOUT_BRAGG_ENERGY              : (  0.00000000,    0.968971,    0.247176)
