"""
Represents a diffraction setup implementation using material data from shadow bragg preprocessor
photon energy in eV
dSpacing returns A
units are in SI.

"""
import numpy

from crystalpy.diffraction.DiffractionSetupAbstract import DiffractionSetupAbstract
from xoppylib.crystals.bragg_preprocessor_file_io import bragg_preprocessor_file_v2_read
from xoppylib.crystals.tools import crystal_fh
# from xoppylib.xoppy_xraylib_util import crystal_fh
# from orangecontrib.xoppy.util.xoppy_xraylib_util import crystal_fh

import scipy.constants as codata

class DiffractionSetupShadowPreprocessorV2(DiffractionSetupAbstract):

    def __init__(self,
                 geometry_type=None,      # Not used, info not in preprocessor file
                 crystal_name="",         # Not used, info not in preprocessor file
                 thickness=1e-6,          # Not used, info not in preprocessor file
                 miller_h=1,              # Not used, info not in preprocessor file
                 miller_k=1,              # Not used, info not in preprocessor file
                 miller_l=1,              # Not used, info not in preprocessor file
                 asymmetry_angle=0.0,     # Not used, info not in preprocessor file
                 azimuthal_angle=0.0,     # Not used, info not in preprocessor file
                 preprocessor_file=""):
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

        self._preprocessor_file = preprocessor_file
        self._preprocessor_dictionary = bragg_preprocessor_file_v2_read(self._preprocessor_file)

    def angleBragg(self, energy):
        """
        Returns the Bragg angle in rad for a given energy.
        :param energy: Energy to calculate the Bragg angle for.
        :return: Bragg angle.
        """
        wavelenth_A = codata.h * codata.c / codata.e / energy * 1e10
        return numpy.arcsin( wavelenth_A / 2 / self.dSpacing())

    def dSpacing(self):
        """
        Returns the lattice spacing d in A.
        :return: Lattice spacing.
        """
        return 1e8 * self._preprocessor_dictionary["dspacing"]

    def unitcellVolume(self):
        """
        Returns the unit cell volume in A**3.

        :return: Unit cell volume
        """
        codata_e2_mc2 = codata.hbar * codata.alpha / codata.m_e / codata.c * 1e2 # in cm
        vol_minusone = self._preprocessor_dictionary["rn"] / codata_e2_mc2
        return 1e24 /vol_minusone


    def F0(self, energy):
        """
        Calculate F0 from Zachariasen.
        :param energy: photon energy in eV.
        :return: F0
        """
        return self.Fall(energy)[0]

    def FH(self, energy, rel_angle=1.0):
        """
        Calculate FH from Zachariasen.
        :param energy: photon energy in eV.
        :return: FH
        """
        return self.Fall(energy, rel_angle=rel_angle)[1]

    def FH_bar(self, energy, rel_angle=1.0):
        """
        Calculate FH_bar from Zachariasen.
        :param energy: photon energy in eV.
        :return: FH_bar
        """
        return self.Fall(energy, rel_angle=rel_angle)[2]

    def Fall(self, energy, rel_angle=1.0):
        # wavelength = codata.h * codata.c / codata.e / energy * 1e10
        # ratio = numpy.sin(self.angleBragg(energy) * rel_angle)/ wavelength
        theta = self.angleBragg(energy) * rel_angle * 1.0
        tmp = crystal_fh(self._preprocessor_dictionary,energy,theta=theta,forceratio=0)

        return tmp['F_0'], tmp['FH'], tmp['FH_BAR']

if __name__ == "__main__":

    from crystalpy.diffraction.GeometryType import BraggDiffraction
    from crystalpy.diffraction.DiffractionSetup import DiffractionSetup
    from xoppylib.crystals.create_bragg_preprocessor_file_v2 import create_bragg_preprocessor_file_v2
    import xraylib
    import numpy


    tmp = create_bragg_preprocessor_file_v2(interactive=False,
                                          DESCRIPTOR="Si", H_MILLER_INDEX=1, K_MILLER_INDEX=1, L_MILLER_INDEX=1,
                                          TEMPERATURE_FACTOR=1.0,
                                          E_MIN=5000.0, E_MAX=15000.0, E_STEP=100.0,
                                          SHADOW_FILE="bragg.dat",
                                          material_constants_library=xraylib)


    a = DiffractionSetupShadowPreprocessorV2(
                 geometry_type=BraggDiffraction,
                 crystal_name="Si", thickness=1e-5,
                 miller_h=1, miller_k=1, miller_l=1,
                 asymmetry_angle=0.0,
                 azimuthal_angle=0.0,
                 preprocessor_file="bragg.dat")

    b = DiffractionSetup(geometry_type=BraggDiffraction,
                 crystal_name="Si", thickness=1e-5,
                 miller_h=1, miller_k=1, miller_l=1,
                 asymmetry_angle=0.0,
                 azimuthal_angle=0.0)

    energy = 8000.0
    print("============ SHADOW / XRAYLIB  ==============")
    print("Photon energy: %g eV " % (energy))
    print("d_spacing: %g %g A " % (a.dSpacing(),b.dSpacing()))
    print("unitCellVolumw: %g %g A**3 " % (a.unitcellVolume(),b.unitcellVolume()))
    print("Bragg angle: %g %g deg " %  (a.angleBragg(energy) * 180 / numpy.pi,
                                     b.angleBragg(energy) * 180 / numpy.pi))
    print("Asymmetry factor b: ", a.asymmetryFactor(energy),
                                b.asymmetryFactor(energy))

    print("F0 ", a.F0(energy), b.F0(energy))
    print("FH ", a.FH(energy), b.FH(energy))
    print("FH_BAR ", a.FH_bar(energy), b.FH_bar(energy))

    print("PSI0 ", a.psi0(energy), b.psi0(energy))
    print("PSIH ", a.psiH(energy), b.psiH(energy))
    print("PSIH_bar ", a.psiH_bar(energy), b.psiH_bar(energy))

    print("DarwinHalfWidths:  ", a.darwinHalfwidth(energy), b.darwinHalfwidth(energy))


    print("\n\n====================== Warning =========================")
    print("Please note a small difference in FH ratio (preprocessor/xraylib): ", a.FH(energy).real /  b.FH(energy).real)
    print("which corresponds to a difference in f0: ")
    print("shadow preprocessor file uses f0_xop() for the coefficients and this is different")
    print("than xraylib.FF_Rayl() by a factor: ")
    ratio = 0.15946847244512372
    import xraylib
    from dabax.dabax_xraylib import DabaxXraylib
    print(DabaxXraylib(file_f0='f0_xop.dat').FF_Rayl(14, 0.15946847244512372) / \
          xraylib.FF_Rayl(14, 0.15946847244512372) )
    print("========================================================\n\n")

    # print("V0: ", a.vectorK0direction(energy).components())
    # print("Bh direction: ", a.vectorHdirection().components())
    # print("Bh: ", a.vectorH().components())
    # print("K0: ", a.vectorK0(energy).components())
    # print("Kh: ", a.vectorKh(energy).components())
    # print("Vh: ", a.vectorKhdirection(energy).components())
    #
    #
    # from crystalpy.util.Photon import Photon
    # print("Difference to ThetaB uncorrected: ",
    #       a.deviationOfIncomingPhoton(Photon(energy_in_ev=energy, direction_vector=a.vectorK0(energy))))
    # #
    # #
    # print("Asymmerey factor b: ", a.asymmetry_factor(energy))
    # print("Bragg angle: %g deg " %  (a.angleBragg(energy) * 180 / numpy.pi))
    # # print("Bragg angle corrected: %g deg " %  (a.angleBraggCorrected(energy) * 180 / numpy.pi))


 #     VIN_BRAGG_UNCORR (Uncorrected): (  0.00000000,    0.968979,   -0.247145)
 #     VIN_BRAGG          (Corrected): (  0.00000000,    0.968971,   -0.247176)
 #     VIN_BRAGG_ENERGY              : (  0.00000000,    0.968971,   -0.247176)
 # Reflected directions matching Bragg angle:
 #    VOUT_BRAGG_UNCORR (Uncorrected): (  0.00000000,    0.968979,    0.247145)
 #    VOUT_BRAGG          (Corrected): (  0.00000000,    0.968971,    0.247176)
 #    VOUT_BRAGG_ENERGY              : (  0.00000000,    0.968971,    0.247176)
