#
# crystal tests
#
import numpy
from dabax.dabax_xraylib import DabaxXraylib

dx = DabaxXraylib(dabax_repository='./dabax/')
print(dx.info())
# dx.set_verbose()

print("DABAX crystal list: \n", dx.Crystal_GetCrystalsList())

crystal = 'Tripalmitin'
siD = dx.Crystal_GetCrystal(crystal)

print("DABAX crystal: \n", siD)

print("%s 111 d-spacing: DABAX: %g" % \
      (crystal, dx.Crystal_dSpacing(siD, 1, 1, 1)))

print("%s 111 bragg angle at 10 keV [deg]: DABAX: %g" % ( crystal, \
    180 / numpy.pi * dx.Bragg_angle(siD, 10, 1, 1, 1), ))

# F0 = dx.Crystal_F_H_StructureFactor(siD,8.0,0,0,0,1.0,ratio_theta_thetaB=1.0)
dabax_all_F = dx.Crystal_F_0_F_H_F_H_bar_StructureFactor(siD, 8.0, 1, 1, 1, 1.0, rel_angle=1.0)

print("F0 dabax: ",
      dx.Crystal_F_H_StructureFactor(siD, 8.0, 0, 0, 0, 1.0, 1.0), dabax_all_F[0],
      )

print("F111 dabax: ",
      dx.Crystal_F_H_StructureFactor(siD, 8.1, 1, 1, 1, 1.0, 1.0), dabax_all_F[1],
      )

print("F-1-1-1 dabax: ",
      dx.Crystal_F_H_StructureFactor(siD, 8.1, -1, -1, -1, 1.0, 1.0), dabax_all_F[2],
      )

print("F001 dabax F2: ", numpy.abs(dx.Crystal_F_H_StructureFactor(siD, 8.1, 0, 0, 1, 1.0, 1.0))**2)
print("F002 dabax F2: ", numpy.abs(dx.Crystal_F_H_StructureFactor(siD, 8.1, 0, 0, 2, 1.0, 1.0))**2)
print("F003 dabax F2: ", numpy.abs(dx.Crystal_F_H_StructureFactor(siD, 8.1, 0, 0, 3, 1.0, 1.0))**2)