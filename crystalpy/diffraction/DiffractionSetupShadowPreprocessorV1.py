"""
Represents a diffraction setup implementation using material data from shadow bragg preprocessor V1
photon energy in eV
dSpacing returns A
units are in SI.
"""
import numpy

from crystalpy.diffraction.DiffractionSetupAbstract import DiffractionSetupAbstract

from crystalpy.util.bragg_preprocessor_file_io import bragg_preprocessor_file_v1_read

import scipy.constants as codata

class DiffractionSetupShadowPreprocessorV1(DiffractionSetupAbstract):
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

    preprocessor_file: str
        The preprocessor file name.

    """

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

        super().__init__(geometry_type=geometry_type,
                         crystal_name=crystal_name,
                         thickness=thickness,
                         miller_h=miller_h,
                         miller_k=miller_k,
                         miller_l=miller_l,
                         asymmetry_angle=asymmetry_angle,
                         azimuthal_angle=azimuthal_angle)

        self._preprocessor_file = preprocessor_file
        self._preprocessor_dictionary = bragg_preprocessor_file_v1_read(self._preprocessor_file)

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
        wavelenth_A = codata.h * codata.c / codata.e / energy * 1e10
        return numpy.arcsin( wavelenth_A / 2 / self.dSpacing())

    def dSpacing(self):
        """Returns the lattice spacing d in A.

        Returns
        -------
        float
            Lattice spacing. in A

        """
        return 1e8 * self._preprocessor_dictionary["dspacing_in_cm"]

    def unitcellVolume(self):
        """Returns the unit cell volume in A^3

        Returns
        -------
        float
            Unit cell volume in A^3.

        """
        codata_e2_mc2 = codata.hbar * codata.alpha / codata.m_e / codata.c * 1e2 # in cm
        vol_minusone = self._preprocessor_dictionary["one_over_volume_times_electron_radius_in_cm"] / codata_e2_mc2
        return 1e24 /vol_minusone


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
        wavelength = codata.h * codata.c / codata.e / energy * 1e10
        ratio = numpy.sin(self.angleBragg(energy) * rel_angle)/ wavelength
        F_0, FH, FH_BAR, STRUCT, FA, FB = self.structure_factor(energy, ratio)
        return F_0, FH, FH_BAR


    #
    #
    #

    def structure_factor(self, energy, ratio):
        """Calculate the structure factors (F_0, FH, FH_BAR, STRUCT, FA, FB)

        Parameters
        ----------
        energy : float or numpy array
            The photon energy in eV.
            
        ratio : float or numpy array
            sin(theta)/lambda

        Returns
        -------
        tuple
            (F_0, FH, FH_BAR, STRUCT, FA, FB)

        """

        F1A, F2A, F1B, F2B = self.__interpolate(energy)
        FOA, FOB = self.__F_elastic(ratio)
        FA = FOA + F1A + 1j * F2A
        FB = FOB + F1B + 1j * F2B

        I_LATT = self._preprocessor_dictionary["i_latt"]

        GA = self._preprocessor_dictionary["ga.real"] + 1j * self._preprocessor_dictionary["ga.imag"]
        GB = self._preprocessor_dictionary["gb.real"] + 1j * self._preprocessor_dictionary["gb.imag"]
        GA_BAR = self._preprocessor_dictionary["ga_bar.real"] + 1j * self._preprocessor_dictionary["ga_bar.imag"]
        GB_BAR = self._preprocessor_dictionary["gb_bar.real"] + 1j * self._preprocessor_dictionary["gb_bar.imag"]
        ATNUM_A = self._preprocessor_dictionary["zeta_a"]
        ATNUM_B = self._preprocessor_dictionary["zeta_b"]
        TEMPER = self._preprocessor_dictionary["temper"]


        # RN = self._preprocessor_dictionary["one_over_volume_times_electron_radius_in_cm"]

        if (I_LATT == 0):
            # ABSORP = 2.0 * RN * R_LAM0 * (4.0*(DIMAG(FA)+DIMAG(FB)))
            F_0 = 4*((F1A + ATNUM_A + F1B + ATNUM_B) + 1j*(F2A + F2B))
        elif (I_LATT == 1):
            # ABSORP = 2.0D0*RN*R_LAM0*(4.0D0*(DIMAG(FA)+DIMAG(FB)))
            F_0 = 4*((F1A + ATNUM_A + F1B + ATNUM_B) + 1j*(F2A + F2B))
        elif (I_LATT == 2):
            FB	 = 0.0 + 0.0j
            # ABSORP = 2.0D0*RN*R_LAM0*(4.0D0*DIMAG(FA))
            F_0 = 4*(F1A + ATNUM_A + 1j*F2A)
        elif (I_LATT == 3):
            # ABSORP = 2.0D0*RN*R_LAM0*(DIMAG(FA)+DIMAG(FB))
            F_0 = (F1A + ATNUM_A + F1B + ATNUM_B) + CI*(F2A + F2B)
        elif (I_LATT == 4):
            FB = 0.0 + 0.0j
            # ABSORP = 2.0D0*RN*R_LAM0*(2.0D0*(DIMAG(FA)))
            F_0 = 2*(F1A+ 1j*F2A)
        elif (I_LATT == 5):
            FB = 0.0 + 0.0j
            # ABSORP = 2.0D0*RN*R_LAM0*(4.0D0*(DIMAG(FA)))
            F_0 = 4*(F1A + 1j*F2A )

        # ! C
        # ! C FH and FH_BAR are the structure factors for (h,k,l) and (-h,-k,-l).
        # ! C
        # ! C srio, Added TEMPER here (95/01/19)
        FH 	= ( (GA * FA) + (GB * FB) ) * TEMPER
        FH_BAR	= ( (GA_BAR * FA) + (GB_BAR * FB) ) * TEMPER
        STRUCT = numpy.sqrt(FH * FH_BAR)
        return F_0, FH, FH_BAR, STRUCT, FA, FB

    def __interpolate(self, PHOT):
        ENERGY = self._preprocessor_dictionary["Energy"]
        NENER = self.__energy_index(PHOT) - 1
        FP_A  = self._preprocessor_dictionary["F1a"]
        FP_B  = self._preprocessor_dictionary["F1b"]
        FPP_A = self._preprocessor_dictionary["F2a"]
        FPP_B = self._preprocessor_dictionary["F2b"]
        F1A	=  FP_A[NENER] +  (FP_A[NENER+1] -  FP_A[NENER]) *  (PHOT - ENERGY[NENER]) / (ENERGY[NENER+1] - ENERGY[NENER])
        F2A	= FPP_A[NENER] + (FPP_A[NENER+1] - FPP_A[NENER]) *  (PHOT - ENERGY[NENER]) / (ENERGY[NENER+1] - ENERGY[NENER])
        F1B	=  FP_B[NENER] +  (FP_B[NENER+1] -  FP_B[NENER]) *  (PHOT - ENERGY[NENER]) / (ENERGY[NENER+1] - ENERGY[NENER])
        F2B	= FPP_B[NENER] + (FPP_B[NENER+1] - FPP_B[NENER]) *  (PHOT - ENERGY[NENER]) / (ENERGY[NENER+1] - ENERGY[NENER])
        return F1A, F2A, F1B, F2B

    def __energy_index(self, energy1):
        Energy = self._preprocessor_dictionary["Energy"]
        energy = numpy.array(energy1)
        if energy.size == 1:
            if (energy < Energy.min()) or (energy > Energy.max()):
                return -100
            ll = numpy.where(Energy > energy)[0][0]
        else:
            ll = numpy.zeros(energy.size, dtype=int)
            for i, ener in enumerate(energy):
                if (ener < Energy.min()) or (ener > Energy.max()):
                    ll[i] = -100
                else:
                    ll[i] = numpy.where( Energy > ener)[0][0]

        NENER = numpy.array(ll)

        if (NENER < 0).any():
            raise Exception("Cannot interpolate: energy outside limits")

        return NENER

    def __F_elastic(self, ratio):
        CA = self._preprocessor_dictionary["fit_a"]
        CB = self._preprocessor_dictionary["fit_b"]
        FOA = CA[2] * ratio ** 2 + CA[1] * ratio + CA[0]
        FOB = CB[2] * ratio ** 2 + CB[1] * ratio + CB[0]
        return FOA, FOB

if __name__ == "__main__":
    if False:
        import numpy
        from crystalpy.diffraction.GeometryType import BraggDiffraction
        from crystalpy.diffraction.DiffractionSetupXraylib import DiffractionSetupXraylib

        try:
            from xoppylib.crystals.create_bragg_preprocessor_file_v1 import create_bragg_preprocessor_file_v1
            import xraylib
            preprocessor_file = "bragg.dat"
            create_bragg_preprocessor_file_v1(interactive=False,
                                                  DESCRIPTOR="Si", H_MILLER_INDEX=1, K_MILLER_INDEX=1, L_MILLER_INDEX=1,
                                                  TEMPERATURE_FACTOR=1.0,
                                                  E_MIN=5000.0, E_MAX=15000.0, E_STEP=100.0,
                                                  SHADOW_FILE=preprocessor_file,
                                                  material_constants_library=xraylib)

        except:
            raise Exception("xoppylib must be installed to create shadow preprocessor files.")



        a = DiffractionSetupShadowPreprocessorV1(
                     geometry_type=BraggDiffraction,
                     crystal_name="Si", thickness=1e-5,
                     miller_h=1, miller_k=1, miller_l=1,
                     asymmetry_angle=0.0,
                     azimuthal_angle=0.0,
                     preprocessor_file=preprocessor_file)

        b = DiffractionSetupXraylib(geometry_type=BraggDiffraction,
                     crystal_name="Si", thickness=1e-5,
                     miller_h=1, miller_k=1, miller_l=1,
                     asymmetry_angle=0.0,
                     azimuthal_angle=0.0)

        energy = 8000.0
        energies = numpy.linspace(energy, energy + 100, 2)

        print("F0 ", a.F0(energy))
        print("F0 [array] ", a.F0(energies))
        print("FH ", a.FH(energy))
        print("FH [array] ", a.FH(energies))
        print("FH_bar ", a.FH_bar(energy))
        print("FH_bar [array] ", a.FH_bar(energies))

        print("============ SHADOW / XRAYLIB  ==============")
        print("Photon energy: %g eV " % (energy))
        print("d_spacing: %g %g A " % (a.dSpacing(),b.dSpacing()))
        print("unitCellVolumw: %g %g A**3 " % (a.unitcellVolume(),b.unitcellVolume()))
        print("Bragg angle: %g %g deg " %  (a.angleBragg(energy) * 180 / numpy.pi,
                                         b.angleBragg(energy) * 180 / numpy.pi))
        print("Asymmetry factor b: ", a.asymmetryFactor(energy),
                                    b.asymmetryFactor(energy))

        print("F0 ", a.F0(energy))
        print("F0 [array] ", a.F0(energies))
        print("FH ", a.FH(energy))
        print("FH [array] ", a.FH(energies))
        print("FH_bar ", a.FH_bar(energy))
        print("FH_bar [array] ", a.FH_bar(energies))

        print("PSI0 ", a.psi0(energy), b.psi0(energy))
        print("PSIH ", a.psiH(energy), b.psiH(energy))
        print("PSIH_bar ", a.psiH_bar(energy), b.psiH_bar(energy))

        print("DarwinHalfWidths:  ", a.darwinHalfwidth(energy), b.darwinHalfwidth(energy))

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

    if True:
        import numpy
        from crystalpy.diffraction.GeometryType import BraggDiffraction
        from crystalpy.diffraction.DiffractionSetupXraylib import DiffractionSetupXraylib

        try:
            from xoppylib.crystals.create_bragg_preprocessor_file_v1 import create_bragg_preprocessor_file_v1
            import xraylib

            tmp = create_bragg_preprocessor_file_v1(interactive=False,
                                                    DESCRIPTOR="Si", H_MILLER_INDEX=1, K_MILLER_INDEX=1,
                                                    L_MILLER_INDEX=1,
                                                    TEMPERATURE_FACTOR=1.0,
                                                    E_MIN=5000.0, E_MAX=15000.0, E_STEP=100.0,
                                                    SHADOW_FILE="bragg.dat",
                                                    material_constants_library=xraylib)
        except:
            raise Exception("xoppylib must be installed to create shadow preprocessor files.")

        a = DiffractionSetupShadowPreprocessorV1(
            geometry_type=BraggDiffraction,
            crystal_name="Si", thickness=1e-5,
            miller_h=1, miller_k=1, miller_l=1,
            asymmetry_angle=0.0,
            azimuthal_angle=0.0,
            preprocessor_file="bragg.dat")

        b = DiffractionSetupXraylib(geometry_type=BraggDiffraction,
                                    crystal_name="Si", thickness=1e-5,
                                    miller_h=1, miller_k=1, miller_l=1,
                                    asymmetry_angle=0.0,
                                    azimuthal_angle=0.0)

        energy = 8000.0
        energies = numpy.linspace(energy, energy + 100, 2)

        print("============ SHADOW / XRAYLIB  ==============")
        print("Photon energy: %g eV " % (energy))
        print("d_spacing: %g %g A " % (a.dSpacing(),
                                       b.dSpacing(),
                                       ))
        print("unitCellVolumw: %g %g A**3 " % (a.unitcellVolume(),
                                               b.unitcellVolume()))
        print("Bragg angle: %g %g deg " % (a.angleBragg(energy) * 180 / numpy.pi,
                                           b.angleBragg(energy) * 180 / numpy.pi))
        print("Bragg angle Corrected: %g %g deg " % (a.angleBraggCorrected(energy) * 180 / numpy.pi,
                                           b.angleBraggCorrected(energy) * 180 / numpy.pi))
        print("Asymmetry factor b: ", a.asymmetryFactor(energy),
              b.asymmetryFactor(energy))

        print("F0 ", a.F0(energy), b.F0(energy))
        print("F0 [array] ", a.F0(energies), b.F0(energies))
        print("FH ", a.FH(energy), b.FH(energy))
        print("FH [array] ", a.FH(energies), b.FH(energies))
        print("FH_bar ", a.FH_bar(energy), b.FH_bar(energy))
        print("FH_bar [array] ", a.FH_bar(energies), b.FH_bar(energies))

        print("PSI0 ", a.psi0(energy), b.psi0(energy))
        print("PSI0  [array] ", a.psi0(energies), b.psi0(energies))
        print("PSIH ", a.psiH(energy), b.psiH(energy))
        print("PSIH  [array] ", a.psiH(energies), b.psiH(energies))
        print("PSIH_bar ", a.psiH_bar(energy), b.psiH_bar(energy))
        print("PSIH_bar  [array] ", a.psiH_bar(energies), b.psiH_bar(energies))

        print("DarwinHalfWidths:  ", a.darwinHalfwidth(energy),
              b.darwinHalfwidth(energy))

        print("\n\n====================== Warning =========================")
        print("Please note a small difference in FH ratio (preprocessor/xraylib): ",
              a.FH(energy).real / b.FH(energy).real)
        print("which corresponds to a difference in f0: ")
        print("shadow preprocessor file uses f0_xop() for the coefficients and this is different")
        print("than xraylib.FF_Rayl() by a factor: ")
        ratio = 0.15946847244512372
        try:
            import xraylib
        except:
            print("xraylib not available")
        from dabax.dabax_xraylib import DabaxXraylib

        print(DabaxXraylib(file_f0='f0_xop.dat').FF_Rayl(14, 0.15946847244512372) / \
              xraylib.FF_Rayl(14, 0.15946847244512372))
        print("========================================================\n\n")

