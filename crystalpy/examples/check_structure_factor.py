from crystalpy.diffraction.GeometryType import BraggDiffraction
from crystalpy.diffraction.DiffractionSetup import DiffractionSetup
from crystalpy.diffraction.DiffractionSetupDabax import DiffractionSetupDabax
from crystalpy.diffraction.DiffractionSetupShadowPreprocessor import DiffractionSetupShadowPreprocessor
from dabax.dabax_xraylib import DabaxXraylib
import numpy

a = DiffractionSetup(geometry_type=BraggDiffraction,
                                       crystal_name="Si", thickness=1e-5,
                                       miller_h=1, miller_k=1, miller_l=1,
                                       asymmetry_angle=0.0,
                                       azimuthal_angle=0.0,)

import socket

if socket.getfqdn().find("esrf") >= 0:
    dx = DabaxXraylib(dabax_repository="http://ftp.esrf.fr/pub/scisoft/DabaxFiles/")
else:
    dx = DabaxXraylib()


a2 = DiffractionSetupDabax(geometry_type=BraggDiffraction,
                                       crystal_name="Si", thickness=1e-5,
                                       miller_h=1, miller_k=1, miller_l=1,
                                       asymmetry_angle=0.0,
                                       azimuthal_angle=0.0,
                                       dabax=dx)


a3 = DiffractionSetupShadowPreprocessor(geometry_type=BraggDiffraction,
                                       crystal_name="Si", thickness=1e-5,
                                       miller_h=1, miller_k=1, miller_l=1,
                                       asymmetry_angle=0.0,
                                       azimuthal_angle=0.0,
                                       preprocessor_file="bragg.dat")

energy = 8000.0
print("Photon energy: %g deg " % (energy))
print("d_spacing: %g %g %g A " %         (a.dSpacing(),a2.dSpacing(),a3.dSpacing()))
print("unitCellVolumw: %g %g %g A**3 " % (a.unitcellVolume(),a2.unitcellVolume(),a3.unitcellVolume()))
print("Bragg angle: %g %g %g deg " % (a.angleBragg(energy) * 180 / numpy.pi,
                                a2.angleBragg(energy) * 180 / numpy.pi,
                                a3.angleBragg(energy) * 180 / numpy.pi))
print("Asymmerey factor b: ", a.asymmetry_factor(energy),
      a2.asymmetry_factor(energy),a3.asymmetry_factor(energy))

print("F0 ",     a.F0(energy), a2.F0(energy), a3.F0(energy))
print("FH ",     a.FH(energy), a2.FH(energy), a3.FH(energy))
print("FH_BAR ", a.FH_bar(energy), a2.FH_bar(energy), a3.FH_bar(energy))

print("PSI0 ",     a.psi0(energy), a2.psi0(energy), a3.psi0(energy))
print("PSIH ",     a.psiH(energy), a2.psiH(energy), a3.psiH(energy))
print("PSIH_bar ", a.psiH_bar(energy), a2.psiH_bar(energy), a3.psiH_bar(energy))