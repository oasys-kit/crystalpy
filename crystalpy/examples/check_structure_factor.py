from crystalpy.diffraction.GeometryType import BraggDiffraction
from crystalpy.diffraction.DiffractionSetupXraylib import DiffractionSetupXraylib
from crystalpy.diffraction.DiffractionSetupDabax import DiffractionSetupDabax
from crystalpy.diffraction.DiffractionSetupShadowPreprocessorV1 import DiffractionSetupShadowPreprocessorV1
from crystalpy.diffraction.DiffractionSetupShadowPreprocessorV2 import DiffractionSetupShadowPreprocessorV2
from dabax.dabax_xraylib import DabaxXraylib
from xoppylib.crystals.create_bragg_preprocessor_file_v1 import create_bragg_preprocessor_file_v1
from xoppylib.crystals.create_bragg_preprocessor_file_v2 import create_bragg_preprocessor_file_v2
import numpy

a = DiffractionSetupXraylib(geometry_type=BraggDiffraction,
                                       crystal_name="Si", thickness=1e-5,
                                       miller_h=1, miller_k=1, miller_l=1,
                                       asymmetry_angle=0.0,
                                       azimuthal_angle=0.0,)


a2 = DiffractionSetupDabax(geometry_type=BraggDiffraction,
                                       crystal_name="Si", thickness=1e-5,
                                       miller_h=1, miller_k=1, miller_l=1,
                                       asymmetry_angle=0.0,
                                       azimuthal_angle=0.0,
                                       dabax=DabaxXraylib())

create_bragg_preprocessor_file_v1(interactive=False,
                                  DESCRIPTOR="Si",
                                  H_MILLER_INDEX=1,
                                  K_MILLER_INDEX=1,
                                  L_MILLER_INDEX=1,
                                  TEMPERATURE_FACTOR=1.0,
                                  E_MIN=5000.0, E_MAX=15000.0, E_STEP=100.0,
                                  SHADOW_FILE="bragg_v1.dat")

a3 = DiffractionSetupShadowPreprocessorV1(geometry_type=BraggDiffraction,
                                       crystal_name="Si", thickness=1e-5,
                                       miller_h=1, miller_k=1, miller_l=1,
                                       asymmetry_angle=0.0,
                                       azimuthal_angle=0.0,
                                       preprocessor_file="bragg_v1.dat")

create_bragg_preprocessor_file_v2(interactive=False,
                                  DESCRIPTOR="Si",
                                  H_MILLER_INDEX=1,
                                  K_MILLER_INDEX=1,
                                  L_MILLER_INDEX=1,
                                  TEMPERATURE_FACTOR=1.0,
                                  E_MIN=5000.0, E_MAX=15000.0, E_STEP=100.0,
                                  SHADOW_FILE="bragg_v2.dat")

a4 = DiffractionSetupShadowPreprocessorV2(geometry_type=BraggDiffraction,
                                       crystal_name="Si", thickness=1e-5,
                                       miller_h=1, miller_k=1, miller_l=1,
                                       asymmetry_angle=0.0,
                                       azimuthal_angle=0.0,
                                       preprocessor_file="bragg_v2.dat")

energy = 8000.0
print("Photon energy: %g deg " % (energy))
print("d_spacing: %g %g %g %g A " %         (a.dSpacing(),a2.dSpacing(),a3.dSpacing(),a4.dSpacing()))
print("unitCellVolumw: %g %g %g %g A**3 " % (a.unitcellVolume(),a2.unitcellVolume(),a3.unitcellVolume(),a4.unitcellVolume()))
print("Bragg angle: %g %g %g %g deg " % (a.angleBragg(energy) * 180 / numpy.pi,
                                a2.angleBragg(energy) * 180 / numpy.pi,
                                a3.angleBragg(energy) * 180 / numpy.pi,
                                a4.angleBragg(energy) * 180 / numpy.pi,
                                         ))

print("F0 ",     a.F0(energy), a2.F0(energy), a3.F0(energy), a4.F0(energy))
print("FH ",     a.FH(energy), a2.FH(energy), a3.FH(energy), a4.FH(energy))
print("FH_BAR ", a.FH_bar(energy), a2.FH_bar(energy), a3.FH_bar(energy), a4.FH_bar(energy))

print("PSI0 ",     a.psi0(energy), a2.psi0(energy), a3.psi0(energy), a4.psi0(energy))
print("PSIH ",     a.psiH(energy), a2.psiH(energy), a3.psiH(energy), a4.psiH(energy))
print("PSIH_bar ", a.psiH_bar(energy), a2.psiH_bar(energy), a3.psiH_bar(energy), a4.psiH_bar(energy))