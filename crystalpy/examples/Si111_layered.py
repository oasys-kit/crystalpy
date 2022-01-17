
#
#
# This example shows the diffraction by a Si 111 crystal calculated in its simplest implementation:
#
#
#    - calculate_simple_diffraction()
#      Uses a crystal setup and calculates the complex transmitivity and reflectivity
#
#
import numpy



from crystalpy.diffraction.GeometryType import BraggDiffraction, BraggTransmission
from crystalpy.diffraction.DiffractionSetup import DiffractionSetup
from crystalpy.diffraction.Diffraction import Diffraction

import scipy.constants as codata

from crystalpy.util.Vector import Vector


from crystalpy.util.Photon import Photon




#
def calculate_simple_diffraction():

    # Create a diffraction setup.

    thickness = 2e-6

    print("\nCreating a diffraction setup...")
    diffraction_setup_r = DiffractionSetup(geometry_type          = BraggDiffraction(),  # GeometryType object
                                               crystal_name           = "Si",                             # string
                                               thickness              = thickness,                             # meters
                                               miller_h               = 1,                                # int
                                               miller_k               = 1,                                # int
                                               miller_l               = 1,                                # int
                                               asymmetry_angle        = 0,#10.0*numpy.pi/180.,            # radians
                                               azimuthal_angle        = 0.0)                              # radians                            # int

    diffraction_setup_t = DiffractionSetup(geometry_type          = BraggTransmission(),  # GeometryType object
                                               crystal_name           = "Si",                             # string
                                               thickness              = thickness,                             # meters
                                               miller_h               = 1,                                # int
                                               miller_k               = 1,                                # int
                                               miller_l               = 1,                                # int
                                               asymmetry_angle        = 0,#10.0*numpy.pi/180.,            # radians
                                               azimuthal_angle        = 0.0)                              # radians

    diffraction_setup_r_half = DiffractionSetup(geometry_type          = BraggDiffraction(),  # GeometryType object
                                               crystal_name           = "Si",                             # string
                                               thickness              = thickness/2,                             # meters
                                               miller_h               = 1,                                # int
                                               miller_k               = 1,                                # int
                                               miller_l               = 1,                                # int
                                               asymmetry_angle        = 0,#10.0*numpy.pi/180.,            # radians
                                               azimuthal_angle        = 0.0)                              # radians                            # int

    diffraction_setup_t_half = DiffractionSetup(geometry_type          = BraggTransmission(),  # GeometryType object
                                               crystal_name           = "Si",                             # string
                                               thickness              = thickness/2,                             # meters
                                               miller_h               = 1,                                # int
                                               miller_k               = 1,                                # int
                                               miller_l               = 1,                                # int
                                               asymmetry_angle        = 0,#10.0*numpy.pi/180.,            # radians
                                               azimuthal_angle        = 0.0)                              # radians



    energy                 = 8000.0                           # eV
    angle_deviation_min    = -300e-6                          # radians
    angle_deviation_max    = 300e-6                           # radians
    angle_deviation_points = 500

    wavelength = codata.h * codata.c / codata.e / energy

    print(">>>>>>>>>>>>", wavelength)
    angle_step = (angle_deviation_max-angle_deviation_min)/angle_deviation_points

    #
    # gets Bragg angle needed to create deviation's scan
    #
    bragg_angle = diffraction_setup_r.angleBragg(energy)

    print("Bragg angle for E=%f eV is %f deg"%(energy,bragg_angle*180.0/numpy.pi))


    # Create a Diffraction object (the calculator)
    diffraction = Diffraction()

    # initialize arrays for storing outputs
    deviations = numpy.zeros(angle_deviation_points)
    intensityS_r = numpy.zeros(angle_deviation_points)
    intensityS_r_half = numpy.zeros(angle_deviation_points)
    intensityS_t = numpy.zeros(angle_deviation_points)

    intensityS_rr = numpy.zeros(angle_deviation_points)
    intensityS_tt = numpy.zeros(angle_deviation_points)

    r = numpy.zeros(angle_deviation_points, dtype=complex)
    r2um = numpy.zeros(angle_deviation_points, dtype=complex)
    t = numpy.zeros(angle_deviation_points, dtype=complex)

    for ia in range(angle_deviation_points):
        deviation = angle_deviation_min + ia * angle_step
        angle = deviation  + bragg_angle

        # calculate the components of the unitary vector of the incident photon scan
        # Note that diffraction plane is YZ
        yy = numpy.cos(angle)
        zz = - numpy.abs(numpy.sin(angle))
        photon = Photon(energy_in_ev=energy,direction_vector=Vector(0.0,yy,zz))

        # perform the calculation
        coeffs_r = diffraction.calculateDiffractedComplexAmplitudes(diffraction_setup_r, photon)
        coeffs_t = diffraction.calculateDiffractedComplexAmplitudes(diffraction_setup_t, photon)
        coeffs_r_half = diffraction.calculateDiffractedComplexAmplitudes(diffraction_setup_r_half, photon)
        coeffs_t_half = diffraction.calculateDiffractedComplexAmplitudes(diffraction_setup_t_half, photon)


        # coeffs_rr = \
        #             coeffs_r_half['S'] * \
        #             (coeffs_t_half['S']**0 + \
        #              coeffs_t_half['S']**2 * ( \
        #                     coeffs_r_half['S'] ** 0 + \
        #                     coeffs_r_half['S'] ** 2 + \
        #                     coeffs_r_half['S'] ** 4 + \
        #                     coeffs_r_half['S'] ** 6 + \
        #                     coeffs_r_half['S'] ** 8 + \
        #                     coeffs_r_half['S'] ** 10 + \
        #                     coeffs_r_half['S'] ** 12 + \
        #                     coeffs_r_half['S'] ** 14 + \
        #                     coeffs_r_half['S'] ** 16 + \
        #                     coeffs_r_half['S'] ** 18 \
        #             ) )

        # a = coeffs_r_half['S']
        # b = coeffs_t_half['S']


        r[ia] = coeffs_r_half['S'].complexAmplitude()
        # t[ia] = coeffs_t_half['S'].complexAmplitude() #* numpy.exp(1j * 2 * numpy.pi / wavelength * (0.5 * thickness / numpy.sin(bragg_angle)) )
        t[ia] = coeffs_t_half['S'].complexAmplitude() * numpy.exp(-1j * 2 * numpy.pi / wavelength * numpy.cos(bragg_angle) * deviation * (thickness/2) )


        r2um[ia] = coeffs_r['S'].complexAmplitude()

        # # sum = a**0
        # # for i in range(2,400,2):
        # #     sum += a**i
        # sum = a**0 / (a**0 - a**2)
        # # coeffs_rr =  a * (b**0 + b**2 * sum)
        # one = a**0
        # coeffs_rr =   a * ( one + b**2 / (one - a**2))
        # coeffs_tt = b**2 * sum
        #
        # intensityS_rr[ia] = coeffs_rr.intensity()
        # intensityS_tt[ia] = coeffs_tt.intensity()
        #
        # # print(coeffs_r)
        # # print(coeffs_r['S'].complexAmplitude())
        #
        # # store results
        deviations[ia] = deviation
        # intensityS_r[ia] = coeffs_r['S'].intensity()
        # intensityS_r_half[ia] = coeffs_r_half['S'].intensity()
        # intensityS_t[ia] = coeffs_t['S'].intensity()


        # print(">>>>>>>>>>", coeffs_r['S'].complexAmplitude() , coeffs_rr.complexAmplitude() )

    # plot results


    print(r, r.shape)

    from srxraylib.plot.gol import plot

    # plot(1e6 * deviations, numpy.abs(r)**2,
    #      1e6 * deviations, numpy.abs(t)**2,
    #     )
    #
    plot(1e6 * deviations, numpy.abs(r)**2,
         1e6 * deviations, numpy.abs(r+r*t*t/(1-r*r))**2,
         1e6 * deviations, numpy.abs(r2um) ** 2,
         1e6 * deviations, numpy.abs(r2um) ** 2 - numpy.abs(r+r*t*t/(1-r*r))**2,
         legend=['r','r2','r 2um','r 2 um - r2']
        )



    # rp = numpy.flip(r)
    # tp = numpy.flip(t)
    # plot(1e6 * deviations, numpy.abs(r)**2,
    #      1e6 * deviations, numpy.abs( r+r*tp*t/(1-r*rp) )**2,
    #      1e6 * deviations, numpy.abs(r2um) ** 2,
    #      legend=['r','r2','r 2um']
    #     )



    # plot(1e6 * deviations, intensityS_r,
    #      1e6 * deviations, intensityS_r_half,
    #      1e6 * deviations, intensityS_rr,
    #      1e6 * deviations, intensityS_tt,
    #      # 1e6 * deviations, intensityS_rr + intensityS_tt,
    #      xtitle="deviation angle [urad]",
    #      ytitle="Reflectivity",
    #      color=['r','b','k'],
    #      linestyle=[None,None,None],
    #      legend=["R 2 um", "R 1 um", "R 1+1um"])

    # plot results
    # from srxraylib.plot.gol import plot
    # plot(1e6*deviations,intensityS_r,
    #      1e6*deviations,intensityS_t,
    #      1e6 * deviations, intensityS_r + intensityS_t,
    #      xtitle="deviation angle [urad]",
    #      ytitle="Reflectivity",
    #      legend=["R Sigma-polarization", "T Sigma-polarization", "R+T"], yrange=[0,1])


#
# main
#
if __name__ == "__main__":


    calculate_simple_diffraction()

