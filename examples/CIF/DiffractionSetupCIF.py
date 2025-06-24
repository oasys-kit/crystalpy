"""
Represents a diffraction setup implementation using crystal data from CIF file and DABAX or XRAYLIB
photon energy in eV
dSpacing returns A
units are in SI.
"""

from crystalpy.diffraction.DiffractionSetupAbstract import DiffractionSetupAbstract
from crystalpy.diffraction.DiffractionSetupDabax import DiffractionSetupDabax
from crystalpy.diffraction.DiffractionSetupXraylib import DiffractionSetupXraylib
from dabax.dabax_xraylib import DabaxXraylib
import xraylib
from CIF2dabax import CIF_GetCrystal

class DiffractionSetupCIFdabax(DiffractionSetupDabax):
    """
    Constructor.

    Parameters
    ----------
    geometry_type: instance of BraggDiffraction, LaueDiffraction, BraggTransmission, or LaueTransmission

    crystal_name: str
        The CIF filename with the crystal.

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

    dabax: DabaxXraylib instance, optional
        Default : None, use DabaxXraylib().

    """

    def __init__(self,
                 geometry_type=None, crystal_name="", thickness=1e-6,
                 miller_h=1, miller_k=1, miller_l=1,
                 asymmetry_angle=0.0,
                 azimuthal_angle=0.0,
                 dabax=None):

        crystal_data = CIF_GetCrystal(crystal_name)
        super().__init__(geometry_type=geometry_type,
                         crystal_name=crystal_name,
                         thickness=thickness,
                         miller_h=miller_h,
                         miller_k=miller_k,
                         miller_l=miller_l,
                         asymmetry_angle=asymmetry_angle,
                         azimuthal_angle=azimuthal_angle,
                         crystal_data=crystal_data,
                         dabax=dabax)

class DiffractionSetupCIFxraylib(DiffractionSetupDabax):
    """
    Constructor.

    Parameters
    ----------
    geometry_type: instance of BraggDiffraction, LaueDiffraction, BraggTransmission, or LaueTransmission

    crystal_name: str
        The CIF filename with the crystal.

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

    """

    def __init__(self,
                 geometry_type=None, crystal_name="", thickness=1e-6,
                 miller_h=1, miller_k=1, miller_l=1,
                 asymmetry_angle=0.0,
                 azimuthal_angle=0.0,
                 ):

        crystal_data = CIF_GetCrystal(crystal_name)
        super().__init__(geometry_type=geometry_type,
                         crystal_name=crystal_name,
                         thickness=thickness,
                         miller_h=miller_h,
                         miller_k=miller_k,
                         miller_l=miller_l,
                         asymmetry_angle=asymmetry_angle,
                         azimuthal_angle=azimuthal_angle,
                         crystal_data=crystal_data)

if __name__ == "__main__":

    if True:
        from crystalpy.diffraction.GeometryType import BraggDiffraction
        import numpy

        # a = DiffractionSetupCIFdabax(geometry_type=BraggDiffraction,
        #         crystal_name = "BiFeO3.cif", thickness = 1e-5,
        #         miller_h = 2, miller_k = 0, miller_l = 4,
        #         asymmetry_angle = 0.0,
        #         azimuthal_angle = 0.0,
        #         dabax = None)

        a = DiffractionSetupCIFxraylib(geometry_type=BraggDiffraction,
                crystal_name="BiFeO3.cif", thickness=1e-5,
                miller_h=2, miller_k=0, miller_l=4,
                asymmetry_angle=0.0,
                azimuthal_angle=0.0)

        energy = 17000.0
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

