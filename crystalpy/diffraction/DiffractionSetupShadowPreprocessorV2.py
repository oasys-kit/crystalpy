"""
Represents a diffraction setup implementation using material data from shadow bragg preprocessor V1
photon energy in eV
dSpacing returns A
units are in SI.
"""
import numpy

from crystalpy.diffraction.DiffractionSetupAbstract import DiffractionSetupAbstract
from crystalpy.util.bragg_preprocessor_file_io import bragg_preprocessor_file_v2_read

import scipy.constants as codata

class DiffractionSetupShadowPreprocessorV2(DiffractionSetupAbstract):
    """
    Constructor.

    Parameters
    ----------
    geometry_type: instance of BraggDiffraction, LaueDiffraction, BraggTransmission, or LaueTransmission
        The crystal geometry.
    crystal_name: str, optional
        The name of the crystal, e.g. "Si".
    thickness: float, optional
        The crystal thickness in m.
    miller_h: int, optional
        Miller index H.
    miller_k: int, optional
        Miller index K.
    miller_l: int, optional
        Miller index L.
    asymmetry_angle: float, optional
        The asymmetry angle between surface normal and Bragg normal (radians).
    azimuthal_angle: float, optional
        The angle between the projection of the Bragg normal on the crystal surface plane and the x axis (radians).
    debye_waller: float, optional
        The Debye-Waller factor exp(-M).
    preprocessor_file: str, optional
        The preprocessor file name.

    Returns
    -------
    instance of DiffractionSetupShadowPreprocessorV2.
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
        self._preprocessor_dictionary = bragg_preprocessor_file_v2_read(self._preprocessor_file)

    def angleBragg(self, energy):
        """
        Returns the Bragg angle for a given energy in radians.

        Parameters
        ----------
        energy : float or numpy array
            Energy to calculate the Bragg angle for. (Default value = 8000.0).

        Returns
        -------
        float or numpy array
            Bragg angle in radians.

        """
        wavelenth_A = codata.h * codata.c / codata.e / energy * 1e10
        return numpy.arcsin( wavelenth_A / 2 / self.dSpacing())

    def dSpacing(self):
        """
        Returns the lattice spacing d in A.

        Returns
        -------
        float
            Lattice spacing in A.
        """
        return 1e8 * self._preprocessor_dictionary["dspacing"]

    def unitcellVolume(self):
        """
        Returns the unit cell volume in A^3

        Returns
        -------
        float
            Unit cell volume in A^3.
        """
        codata_e2_mc2 = codata.hbar * codata.alpha / codata.m_e / codata.c * 1e2 # in cm
        vol_minusone = self._preprocessor_dictionary["rn"] / codata_e2_mc2
        return 1e24 /vol_minusone


    def F0(self, energy):
        """
        Calculate the structure factor F0.

        Parameters
        ----------
        energy : float or numpy array
            photon energy in eV. (Default value = 8000.0).

        Returns
        -------
        complex or numpy array
            F0.
        """
        return self.Fall(energy)[0]

    def FH(self, energy, rel_angle=1.0):
        """
        Calculate the structure factor FH.

        Parameters
        ----------
        energy : float or numpy array
            photon energy in eV.
        rel_angle: float or numpy array, optional
            ratio of the incident angle and the Bragg angle (Default : rel_angle=1.0)

        Returns
        -------
        complex or numpy array
            FH.
        """
        return self.Fall(energy, rel_angle=rel_angle)[1]

    def FH_bar(self, energy, rel_angle=1.0):
        """
        Calculate the structure factor  FH_bar.

        Parameters
        ----------
        energy : float or numpy array
            photon energy in eV.
        rel_angle: float or numpy array, optional
            ratio of the incident angle and the Bragg angle (Default : rel_angle=1.0).

        Returns
        -------
        complex or numpy array
            FH_bar.
        """
        return self.Fall(energy, rel_angle=rel_angle)[2]

    def Fall(self, energy, rel_angle=1.0):
        """
        Calculate the all structure factor  (F0, FH, FH_bar).

        Parameters
        ----------
        energy : float or numpy array
            photon energy in eV.
        rel_angle: float or numpy array, optional
            ratio of the incident angle and the Bragg angle (Default : rel_angle=1.0).

        Returns
        -------
        tuple
            (F0, FH, FH_bar).
        """
        # wavelength = codata.h * codata.c / codata.e / energy * 1e10
        # ratio = numpy.sin(self.angleBragg(energy) * rel_angle)/ wavelength
        theta = self.angleBragg(energy) * rel_angle * 1.0
        tmp = crystal_fh(self._preprocessor_dictionary, energy, theta=theta, forceratio=0)
        return tmp['F_0'], tmp['FH'], tmp['FH_BAR']


# this function is copied from
# https://github.com/oasys-kit/xoppylib/blob/main/xoppylib/crystals/tools.py
# to make crystalpy not depending on xoppylib
# Vectorized by SSLS: Xiaojiang Yu, xiaojiang@nus.edu.sg
# todo: copy tthis vectorized version back to xoppylib
def crystal_fh(input_dictionary, phot_in, theta=None, forceratio=0):
    """
    Computes the structure factors and other parameters for a given photon energy (or array).

    Parameters
    ----------
    input_dictionary : dict
        as resulting from bragg_calc().
    phot_in : float or numpy array
        photon energy in eV.
    theta : float or numpy array, optional
        incident angle (half of scattering angle) in rad (Default value = None).
    forceratio : int, float or numpy array, optional
         Value of sin(theta)/lambda to be used when calculating F0. By default, forceratio=0, it means
         that this ratio is automatically calculated used the Bragg angle.

    Returns
    -------
    dict
        return {"PHOT":phot, "WAVELENGTH":r_lam0*1e-2 ,"THETA":itheta, "F_0":F_0, "FH":FH, "FH_BAR":FH_BAR,
        "STRUCT":STRUCT, "psi_0":psi_0, "psi_h":psi_h, "psi_hbar":psi_hbar,
        "DELTA_REF":DELTA_REF, "REFRAC":REFRAC, "ABSORP":ABSORP, "RATIO":ratio,
        "ssr":ssr, "spr":spr, "psi_over_f":psi_over_f, "info":txt}.
    """

    # outfil    = input_dictionary["outfil"]
    # fract     = input_dictionary["fract"]
    rn        = input_dictionary["rn"]
    dspacing  = numpy.array(input_dictionary["dspacing"])
    nbatom    = numpy.array(input_dictionary["nbatom"])
    atnum     = numpy.array(input_dictionary["atnum"])
    temper    = numpy.array(input_dictionary["temper"])
    G_0       = numpy.array(input_dictionary["G_0"])
    G         = numpy.array(input_dictionary["G"])
    G_BAR     = numpy.array(input_dictionary["G_BAR"])
    f0coeff   = numpy.array(input_dictionary["f0coeff"])
    npoint    = numpy.array(input_dictionary["npoint"])
    energy    = numpy.array(input_dictionary["energy"])
    fp        = numpy.array(input_dictionary["f1"])
    fpp       = numpy.array(input_dictionary["f2"])
    fraction  = numpy.array(input_dictionary["fraction"])


    phot_in = numpy.array(phot_in, dtype=float).reshape(-1)
    theta   = numpy.array(theta,   dtype=float).reshape(-1)

    toangstroms = codata.h * codata.c / codata.e * 1e10

    #To fully exploit numpy's powerful calculation of matrix
    #SSLS: Xiaojiang Yu, xiaojiang@nus.edu.sg
    if ( (phot_in < energy[0]).any() or (phot_in > energy[-1]).any()):
            raise Exception("Photon energy outside of valid limits [%g,%g]" % (energy[0], energy[-1]))

    if theta is None:
        itheta = numpy.arcsin(toangstroms * 1e-8 / phot_in / 2 / dspacing)
    else:
        itheta = theta

    if forceratio == 0:
        ratio = numpy.sin(itheta) / (toangstroms / phot_in)
    else:
        ratio = numpy.array([1 / (2 * dspacing * 1e8)] * len(phot_in))
    # print("Ratio: ",ratio)

    F0   = numpy.zeros((len(phot_in), nbatom))
    F000 = numpy.zeros((len(phot_in), nbatom))
    for j in range(nbatom):
        #icentral = int(f0coeff.shape[1]/2)
        #F0[j] = f0coeff[j,icentral]
        icentral = int(len(f0coeff[j]) / 2)
        F0[:,j] = f0coeff[j][icentral]
        # F000[j] = F0[j]
        for i in range(icentral):
            #F0[j] += f0coeff[j,i] * numpy.exp(-1.0*f0coeff[j,i+icentral+1]*ratio**2)
            F0[:,j] += f0coeff[j][i] * numpy.exp(-1.0 * f0coeff[j][i + icentral + 1] * ratio**2)
            #srio F000[j] += f0coeff[j][i]  #actual number of electrons carried by each atom, X.J. Yu, slsyxj@nus.edu.sg
        F000[:, j] = atnum[j] # srio
    # ;C
    # ;C Interpolate for the atomic scattering factor.
    # ;C

    nener = numpy.zeros((len(phot_in), 2), dtype=int)
    L = len(energy) - 1
    for j, ienergy in enumerate(phot_in):
        test = energy <= ienergy
        nener[j, 0] = numpy.where(test==True)[0][-1]
        if (nener[j, 0] >= L):
            nener[j, 1] = L
        else:
            nener[j,1] = nener[j, 0] + 1

    F1 = numpy.zeros((len(phot_in),nbatom),dtype=float)
    F2 = numpy.zeros((len(phot_in),nbatom),dtype=float)
    F  = numpy.zeros((len(phot_in),nbatom),dtype=complex)

    for j in range(nbatom):
        F1[:,j] = fp[j,nener[:, 0]] + (fp[j,nener[:, 1]] - fp[j,nener[:, 0]]) * \
        (phot_in - energy[nener[:, 0]]) / (energy[nener[:, 1]] - energy[nener[:, 0]])
        F2[:,j] = fpp[j,nener[:, 0]] + (fpp[j,nener[:, 1]] - fpp[j,nener[:, 0]]) * \
        (phot_in - energy[nener[:, 0]]) / (energy[nener[:, 1]] - energy[nener[:, 0]])

    r_lam0 = toangstroms * 1e-8 / phot_in
    for j in range(nbatom):
        F[:,j] = F0[:,j] + F1[:,j] + 1j * F2[:,j]
        # print("F",F)

    F_0 = numpy.linspace(0, 0, len(phot_in), dtype=complex)
    FH      = F_0.copy()
    FH_BAR  = F_0.copy()
    FHr     = F_0.copy()
    FHi     = F_0.copy()
    FH_BARr = F_0.copy()
    FH_BARi = F_0.copy()

    TEMPER_AVE = 1.0
    for j in range(nbatom):
        FH  += fraction[j] * (G[j] *   F[:, j] * 1.0) * temper[j]
        FHr += fraction[j] * (G[j] *  (F0[:, j] + F1[:, j])* 1.0) * temper[j]
        FHi += fraction[j] * (G[j] *   F2[:, j] * 1.0) * temper[j]
        FN = F000[:, j] + F1[:, j] + 1j * F2[:, j]
        F_0 += fraction[j] * (G_0[j] *  FN  * 1.0)
        # TEMPER_AVE *= (temper[j])**(G_0[j]/(G_0.sum()))

        FH_BAR  += fraction[j] * ((G_BAR[j] * F[:, j] * 1.0)) * temper[j]
        FH_BARr += fraction[j] * ((G_BAR[j] * (F0[:, j]  + F1[:, j]) * 1.0)) * temper[j]
        FH_BARi += fraction[j] * ((G_BAR[j] * F2[:, j] * 1.0)) * temper[j]
        # print("TEMPER_AVE: ",TEMPER_AVE)

    # ;C
    # ;C multiply by the average temperature factor
    # ;C

    # FH      *= TEMPER_AVE
    # FHr     *= TEMPER_AVE
    # FHi     *= TEMPER_AVE
    # FH_BAR  *= TEMPER_AVE
    # FH_BARr *= TEMPER_AVE
    # FH_BARi *= TEMPER_AVE

    STRUCT = numpy.sqrt(FH * FH_BAR)

    # ;C
    # ;C   PSI_CONJ = F*( note: PSI_HBAR is PSI at -H position and is
    # ;C   proportional to fh_bar but PSI_CONJ is complex conjugate os PSI_H)
    # ;C


    psi_over_f = rn * r_lam0**2 / numpy.pi
    psi_h      = rn * r_lam0**2 / numpy.pi * FH
    psi_hr     = rn * r_lam0**2 / numpy.pi * FHr
    psi_hi     = rn * r_lam0**2 / numpy.pi * FHi
    psi_hbar   = rn * r_lam0**2 / numpy.pi * FH_BAR
    psi_hbarr  = rn * r_lam0**2 / numpy.pi * FH_BARr
    psi_hbari  = rn * r_lam0**2 / numpy.pi * FH_BARi
    psi_0      = rn * r_lam0**2 / numpy.pi * F_0
    psi_conj   = rn * r_lam0**2 / numpy.pi * FH.conjugate()

    # ;
    # ; Darwin width
    # ;
    # print(rn,r_lam0,STRUCT,itheta)
    ssvar = rn * (r_lam0**2) * STRUCT / numpy.pi / numpy.sin(2.0 * itheta)
    spvar = ssvar * numpy.abs((numpy.cos(2.0 * itheta)))
    ssr = ssvar.real
    spr = spvar.real

    # ;C
    # ;C computes refractive index.
    # ;C ([3.171] of Zachariasen's book)
    # ;C
    REFRAC = (1.0+0j) - r_lam0**2 * rn * F_0 / 2/ numpy.pi
    DELTA_REF = 1.0 - REFRAC.real
    ABSORP = 4.0 * numpy.pi * (-REFRAC.imag) / r_lam0

    #only output the first calculation value, however, other value is possible
    #otherwise will create huge txt content for large number of ray number
    #SSLS: Xiaojiang Yu
    txt = ""
    txt += '\n******************************************************'
    txt += '\n       at energy    = '+repr(phot_in[0])+' eV'
    txt += '\n                    = '+repr(r_lam0[0]*1e8)+' Angstroms'
    txt += '\n       and at angle = '+repr(itheta[0]*180.0 / numpy.pi) + ' degrees'
    txt += '\n                    = '+repr(itheta[0])+' rads'
    txt += '\n******************************************************'

    for j in range(nbatom):
        txt += '\n  '
        txt += '\nFor atom '+repr(j+1)+':'
        txt += '\n       fo + fp+ i fpp = '
        txt += '\n        '+repr(F0[0,j])+' + '+ repr(F1[0,j].real)+' + i'+ repr(F2[0,j])+" ="
        txt += '\n        '+repr(F0[0,j] + F1[0,j] + 1j * F2[0,j])
        txt += '\n       Z = '+repr(atnum[j])
        txt += '\n       Temperature factor = '+repr(temper[j])
    txt += '\n  '
    txt += '\n Structure factor F(0,0,0) = '+repr(F_0[0])
    txt += '\n Structure factor FH = '      +repr(FH[0])
    txt += '\n Structure factor FH_BAR = '  +repr(FH_BAR[0])
    txt += '\n Structure factor F(h,k,l) = '+repr(STRUCT[0])
    txt += '\n  '
    txt += '\n Psi_0  = '   +repr(psi_0[0])
    txt += '\n Psi_H  = '   +repr(psi_h[0])
    txt += '\n Psi_HBar  = '+repr(psi_hbar[0])
    txt += '\n  '
    txt += '\n Psi_H(real) Real and Imaginary parts = '   + repr(psi_hr[0])
    txt += '\n Psi_H(real) Modulus  = '                   + repr(numpy.abs(psi_hr[0]))
    txt += '\n Psi_H(imag) Real and Imaginary parts = '   + repr(psi_hi[0])
    txt += '\n Psi_H(imag) Modulus  = '                   + repr(abs(psi_hi[0]))
    txt += '\n Psi_HBar(real) Real and Imaginary parts = '+ repr(psi_hbarr[0])
    txt += '\n Psi_HBar(real) Modulus  = '                + repr(abs(psi_hbarr[0]))
    txt += '\n Psi_HBar(imag) Real and Imaginary parts = '+ repr(psi_hbari[0])
    txt += '\n Psi_HBar(imag) Modulus  = '                + repr(abs(psi_hbari[0]))
    txt += '\n  '
    txt += '\n Psi/F factor = '                           + repr(psi_over_f[0])
    txt += '\n  '
    txt += '\n Average Temperature factor = '             + repr(TEMPER_AVE)
    txt += '\n Refraction index = 1 - delta - i*beta'
    txt += '\n            delta = '                       + repr(DELTA_REF[0])
    txt += '\n             beta = '                       + repr(1.0e0*REFRAC[0].imag)
    txt += '\n Absorption coeff = '                       + repr(ABSORP[0])+' cm^-1'
    txt += '\n  '
    txt += '\n e^2/(mc^2)/V = '                           + repr(rn)+' cm^-2'
    txt += '\n d-spacing = '                              + repr(dspacing*1.0e8)+' Angstroms'
    txt += '\n SIN(theta)/Lambda = '                      + repr(ratio[0])
    txt += '\n  '
    txt += '\n Darwin width for symmetric s-pol [microrad] = ' + repr(2.0e6*ssr[0])
    txt += '\n Darwin width for symmetric p-pol [microrad] = ' + repr(2.0e6*spr[0])

    return {"PHOT":phot_in, "WAVELENGTH":r_lam0*1e-2 ,"THETA":itheta, "F_0":F_0, "FH":FH, "FH_BAR":FH_BAR,
	        "STRUCT":STRUCT, "psi_0":psi_0, "psi_h":psi_h, "psi_hbar":psi_hbar,
        	"DELTA_REF":DELTA_REF, "REFRAC":REFRAC, "ABSORP":ABSORP, "RATIO":ratio,
        	"ssr":ssr, "spr":spr, "psi_over_f":psi_over_f, "info":txt}

if __name__ == "__main__":
    if 1:
        import numpy
        from crystalpy.diffraction.GeometryType import BraggDiffraction
        from crystalpy.diffraction.DiffractionSetupXraylib import DiffractionSetupXraylib

        try:
            from xoppylib.crystals.create_bragg_preprocessor_file_v2 import create_bragg_preprocessor_file_v2
            import xraylib
            tmp = create_bragg_preprocessor_file_v2(interactive=False,
                                                  DESCRIPTOR="Si", H_MILLER_INDEX=1, K_MILLER_INDEX=1, L_MILLER_INDEX=1,
                                                  TEMPERATURE_FACTOR=1.0,
                                                  E_MIN=5000.0, E_MAX=15000.0, E_STEP=100.0,
                                                  SHADOW_FILE="bragg.dat",
                                                  material_constants_library=xraylib)
        except:
            raise Exception("xoppylib must be installed to create shadow preprocessor files.")


        a = DiffractionSetupShadowPreprocessorV2(
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
        print("Bragg angle: %g %g deg " %  (a.angleBragg(energy) * 180 / numpy.pi,
                                            b.angleBragg(energy) * 180 / numpy.pi))
        print("Asymmetry factor b: ", a.asymmetryFactor(energy),
                                      b.asymmetryFactor(energy))

        print("F0 ",             a.F0(energy),       b.F0(energy))
        print("F0 [array] ",     a.F0(energies),     b.F0(energies))
        print("FH ",             a.FH(energy),       b.FH(energy))
        print("FH [array] ",     a.FH(energies),     b.FH(energies))
        print("FH_bar ",         a.FH_bar(energy),   b.FH_bar(energy))
        print("FH_bar [array] ", a.FH_bar(energies), b.FH_bar(energies))

        print("PSI0 ",     a.psi0(energy),     b.psi0(energy))
        print("PSIH ",     a.psiH(energy),     b.psiH(energy))
        print("PSIH_bar ", a.psiH_bar(energy), b.psiH_bar(energy))

        print("DarwinHalfWidths:  ", a.darwinHalfwidth(energy),
                                     b.darwinHalfwidth(energy))


        print("\n\n====================== Warning =========================")
        print("Please note a small difference in FH ratio (preprocessor/xraylib): ", a.FH(energy).real /  b.FH(energy).real)
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
        # print("Asymmerey factor b: ", a.asymmetryFactor(energy))
        # print("Bragg angle: %g deg " %  (a.angleBragg(energy) * 180 / numpy.pi))
        # print("Bragg angle corrected: %g deg " %  (a.angleBraggCorrected(energy) * 180 / numpy.pi))
        #
        #
        #  #     VIN_BRAGG_UNCORR (Uncorrected): (  0.00000000,    0.968979,   -0.247145)
        #  #     VIN_BRAGG          (Corrected): (  0.00000000,    0.968971,   -0.247176)
        #  #     VIN_BRAGG_ENERGY              : (  0.00000000,    0.968971,   -0.247176)
        #  # Reflected directions matching Bragg angle:
        #  #    VOUT_BRAGG_UNCORR (Uncorrected): (  0.00000000,    0.968979,    0.247145)
        #  #    VOUT_BRAGG          (Corrected): (  0.00000000,    0.968971,    0.247176)
        #  #    VOUT_BRAGG_ENERGY              : (  0.00000000,    0.968971,    0.247176)