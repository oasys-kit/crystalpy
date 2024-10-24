
import numpy



from crystalpy.diffraction.GeometryType import LaueDiffraction, LaueTransmission, BraggDiffraction, BraggTransmission
from crystalpy.diffraction.DiffractionSetupXraylib import DiffractionSetupXraylib
from crystalpy.diffraction.Diffraction import Diffraction


from crystalpy.util.Vector import Vector


from crystalpy.util.Photon import Photon


from srxraylib.plot.gol import plot_image, plot


#
def calculate_diffraction_map_inside_crystal(geometry_type=LaueDiffraction(),
                                             asymmetry_angle=numpy.radians(65),
                                             thickness=10e-6,
                                             calculation_strategy_flag=0):

    # Create a diffraction setup.

    print("\nCreating a diffraction setup...")


    diffraction_setup = DiffractionSetupXraylib(geometry_type=geometry_type,  # GeometryType object
                                         crystal_name="Si",  # string
                                         thickness=thickness,  # meters
                                         miller_h=1,  # int
                                         miller_k=1,  # int
                                         miller_l=1,  # int
                                         asymmetry_angle= asymmetry_angle,  # 10.0*numpy.pi/180.,            # radians
                                         azimuthal_angle=0.0)  # radians                            # int
    # int


    energy                 = 8000.0                           # eV
    angle_deviation_min    = -100e-6                          # radians
    angle_deviation_max    =  100e-6                           # radians
    angle_deviation_points = 100

    s_ratio_deviation_points = 50

    angle_step = (angle_deviation_max-angle_deviation_min)/angle_deviation_points

    #
    # gets Bragg angle needed to create deviation's scan
    #
    bragg_angle = diffraction_setup.angleBragg(energy)

    print("Bragg angle for E=%f eV is %f deg"%(energy,bragg_angle*180.0/numpy.pi))


    # initialize arrays for storing outputs
    deviations = numpy.zeros(angle_deviation_points)


    s_ratios = numpy.linspace(-1.0,0.0,s_ratio_deviation_points)

    import time

    t0 = time.time()

    YY = numpy.zeros( angle_deviation_points *  s_ratio_deviation_points, dtype=float)
    ZZ = numpy.zeros_like(YY)
    S_RATIOS = numpy.zeros_like(YY)

    # iii = -1
    # for ia in range(angle_deviation_points):
    #     for isr in range(s_ratio_deviation_points):
    #         s_ratio = -s_ratios[isr]
    #         iii += 1
    #         deviation = angle_deviation_min + ia * angle_step
    #         angle = deviation  + bragg_angle + asymmetry_angle
    #         # calculate the components of the unitary vector of the incident photon scan
    #         # Note that diffraction plane is YZ
    #         yy = numpy.cos(angle)
    #         zz = - numpy.abs(numpy.sin(angle))
    #         YY[iii] = yy
    #         ZZ[iii] = zz
    #         S_RATIOS[iii] = s_ratio
    #         photon = Photon(energy_in_ev=energy,direction_vector=Vector(0.0,yy,zz))
    #
    #         if False:
    #             # perform the calculation
    #             coeffs = Diffraction.calculateDiffractedComplexAmplitudes(diffraction_setup, photon,
    #                                                                       calculation_method=1,
    #                                                                       calculation_strategy_flag=calculation_strategy_flag)
    #         else:
    #             # Get PerfectCrystal instance for the current photon.
    #             perfect_crystal = Diffraction._perfectCrystalForPhoton(diffraction_setup, photon,
    #                                                                    calculation_strategy_flag=calculation_strategy_flag)
    #             coeffs = perfect_crystal.calculateDiffractionGuigay(photon, debug=0, s_ratio=s_ratio)
    #
    #         # store results
    #         deviations[ia] = deviation
    #         intensityS[ia, isr] = numpy.abs(coeffs['S']) ** 2
    #         intensityP[ia, isr] = numpy.abs(coeffs['P']) ** 2


    t1 = time.time()
    print('>> full loop: ', t1 - t0 )

    iii = -1
    for ia in range(angle_deviation_points):
        for isr in range(s_ratio_deviation_points):
            s_ratio = -s_ratios[isr]
            iii += 1
            deviation = angle_deviation_min + ia * angle_step
            angle = deviation  + bragg_angle + asymmetry_angle
            # calculate the components of the unitary vector of the incident photon scan
            # Note that diffraction plane is YZ
            yy = numpy.cos(angle)
            zz = - numpy.abs(numpy.sin(angle))
            YY[iii] = yy
            ZZ[iii] = zz
            S_RATIOS[iii] = s_ratio
            deviations[ia] = deviation

    t2 = time.time()
    print('>> feeding loop: ', t2 - t1 )


    PHOTON = Photon(energy_in_ev=numpy.ones_like(YY) * energy, direction_vector=Vector(0.0*YY,YY,ZZ))

    Perfect_crystal = Diffraction._perfectCrystalForPhoton(diffraction_setup, PHOTON,
                                                           calculation_strategy_flag=calculation_strategy_flag)
    Coeffs = Perfect_crystal.calculateDiffractionGuigay(PHOTON, debug=0, s_ratio=S_RATIOS)

    t3 = time.time()
    print('>> vector loop: ', t3 - t2, t3 - t1 )

    print(Coeffs['S'].shape, Coeffs['S'][0])

    INTENSITIES =  numpy.abs(numpy.array(Coeffs['S'], dtype=complex)) ** 2
    INTENSITIES =  numpy.reshape(INTENSITIES, (angle_deviation_points, s_ratio_deviation_points))
    # s_ext = 0.5 / (numpy.abs(Coeffs['a'].imag - 2 * Coeffs['u0'].imag))
    s_ext = 0.5 / (numpy.abs(Coeffs['a'].imag ))
    S_EXT = numpy.reshape(s_ext, (angle_deviation_points, s_ratio_deviation_points))
    #
    return INTENSITIES, deviations, s_ratios, S_EXT


#
# main
#
if __name__ == "__main__":
    calculation_strategy_flag = 0  # 0=mpmath 1=numpy 2=numpy-truncated
    geometry_type = BraggTransmission() # BraggDiffraction()
    thickness = 50e-6
    intensityS, deviations, s_ratios, S_EXT = calculate_diffraction_map_inside_crystal(geometry_type=geometry_type, asymmetry_angle=numpy.radians(0), thickness=thickness, calculation_strategy_flag=calculation_strategy_flag)


    if geometry_type == LaueDiffraction():
        title = r"Laue $|D_H|^2$"
        iii = 0
        filename='Laue_DH.png'
    elif geometry_type == LaueTransmission():
        title = r"Laue $|D_0|^2$"
        iii = 0
        filename='Laue_D0.png'
    elif geometry_type == BraggDiffraction():
        title = r"Bragg $|D_H|^2$"
        iii = -1
        filename='Bragg_DH.png'
    elif geometry_type == BraggTransmission():
        title = r"Bragg $|D_0|^2$"
        iii = 0
        filename='Bragg_D0.png'

    plot_image(intensityS, deviations*1e6, s_ratios, xtitle=r"$\theta-\theta_B$ [$\mu$rad]", ytitle="-s/T", title=title, aspect='auto', show=0)

    plot_image(intensityS, deviations * 1e6, -s_ratios * thickness * 1e6, xtitle=r"$\theta-\theta_B$ [$\mu$m]", ytitle="depth [um]", title=title,
               aspect='auto', show=0)

    plot(deviations*1e6, intensityS[:,iii], title=iii, xtitle=r"$\theta-\theta_B$ [$\mu$rad]", ytitle="I_0 at back surface", show=0)


    nMU = deviations.size
    EXT = numpy.zeros(nMU)
    EXT_THEORY = numpy.zeros(nMU)

    depth_in_um = -s_ratios * thickness * 1e6

    print(">>>> S_EXT", S_EXT.shape )
    for i in range(nMU):
        # plot(depth_in_um, intensityS[i, :])
        ylog = numpy.log( intensityS[i, :] )
        # P = numpy.polyfit(depth_in_um, numpy.log(intensityS[i, :]), 1 )
        P = numpy.polyfit(depth_in_um, numpy.log(intensityS[i, :]), 1 )
        EXT[i] = -1.0 / P[0]
        EXT_THEORY[i] = S_EXT[i, 0] * 1e6 * numpy.sin(numpy.radians(14.308608)) # intensity
        print(">>>>", i, S_EXT[i, :])
        # plot(depth_in_um, intensityS[i, :],
        #      depth_in_um, numpy.exp(-depth_in_um/ (EXT[i])),
        #      legend=['data','fit %d' % i])


    import matplotlib
    import matplotlib.pylab as plt
    matplotlib.rcParams.update({'font.size': 14})
    plot(
         deviations * 1e6, EXT,
        deviations * 1e6, EXT_THEORY,
         xtitle=r"$\theta-\theta_B$ [$\mu$rad]", ytitle=r"Extinction depth [$\mu$m]",
         legend=['Fitted', 'Calculated', ],
         marker  =['+',None,],
         linestyle=['', None,],
         title="",
         show=0)

    plt.savefig('penetration.pdf')
    plt.show()