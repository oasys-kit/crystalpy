
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



from crystalpy.diffraction.GeometryType import BraggDiffraction
from crystalpy.diffraction.DiffractionSetup import DiffractionSetup
from crystalpy.diffraction.DiffractionSetupDabax import DiffractionSetupDabax
from crystalpy.diffraction.Diffraction1 import Diffraction1 as Diffraction


from crystalpy.util.Vector import Vector
from crystalpy.util.Photon import Photon
from crystalpy.diffraction.PerfectCrystalDiffraction import PerfectCrystalDiffraction


from dabax.dabax_xraylib import DabaxXraylib

import time


def calculate_simple_diffraction_angular_scan_accelerated(calculation_method=0):

    # Create a diffraction setup.


    diffraction_setup_dabax = DiffractionSetupDabax(geometry_type          = BraggDiffraction(),  # GeometryType object
                                               crystal_name           = "YB66",                             # string
                                               thickness              = 7e-3,                             # meters
                                               miller_h               = 4,                                # int
                                               miller_k               = 0,                                # int
                                               miller_l               = 0,                                # int
                                               asymmetry_angle        = 0,#10.0*numpy.pi/180.,                              # radians
                                               azimuthal_angle        = 0.0,
                                               dabax=DabaxXraylib())                              # radians


    energy                 = 8040.0                           # eV
    angle_deviation_min    = 20e-6                          # radians
    angle_deviation_max    = 80e-6                           # radians
    angle_deviation_points = 200

    angle_step = (angle_deviation_max-angle_deviation_min)/angle_deviation_points

    #
    # gets Bragg angle needed to create deviation's scan
    #
    bragg_angle = diffraction_setup_dabax.angleBragg(energy)

    print("Bragg angle for E=%f eV is %f deg"%(energy,bragg_angle*180.0/numpy.pi))


    # Create a Diffraction object (the calculator)
    diffraction = Diffraction()
    diffraction_dabax = Diffraction()

    # initialize arrays for storing outputs
    deviations = numpy.zeros(angle_deviation_points)
    intensityS = numpy.zeros(angle_deviation_points)
    intensityP = numpy.zeros(angle_deviation_points)
    intensityS_dabax = numpy.zeros(angle_deviation_points)
    intensityP_dabax = numpy.zeros(angle_deviation_points)

    for ia in range(angle_deviation_points):
        deviation = angle_deviation_min + ia * angle_step
        angle = deviation  + bragg_angle

        # calculate the components of the unitary vector of the incident photon scan
        # Note that diffraction plane is YZ
        yy = numpy.cos(angle)
        zz = - numpy.abs(numpy.sin(angle))
        photon = Photon(energy_in_ev=energy,direction_vector=Vector(0.0,yy,zz))

        # perform the calculation
        coeffs = diffraction.calculateDiffractedComplexAmplitudes(diffraction_setup_dabax, photon, method=calculation_method)

        # store results
        deviations[ia] = deviation
        intensityS[ia] = coeffs['S'].intensity()
        intensityP[ia] = coeffs['P'].intensity()

    psi_0, psi_H, psi_H_bar = diffraction_setup_dabax.psiAll(energy)

    for ia in range(angle_deviation_points):
        deviation = angle_deviation_min + ia * angle_step
        angle = deviation  + bragg_angle

        # calculate the components of the unitary vector of the incident photon scan
        # Note that diffraction plane is YZ
        yy = numpy.cos(angle)
        zz = - numpy.abs(numpy.sin(angle))
        photon = Photon(energy_in_ev=energy,direction_vector=Vector(0.0,yy,zz))

        # perform the calculation
        # coeffs_dabax = diffraction_dabax.calculateDiffractedComplexAmplitudes(diffraction_setup_dabax, photon)
        #
        # # store results
        # deviations[ia] = deviation
        # intensityS_dabax[ia] = coeffs_dabax['S'].intensity()
        # intensityP_dabax[ia] = coeffs_dabax['P'].intensity()


        # Create PerfectCrystalDiffraction instance.
        perfect_crystal = PerfectCrystalDiffraction(geometry_type=diffraction_setup_dabax.geometryType(),
                                                    bragg_normal=diffraction_setup_dabax.normalBragg(),
                                                    surface_normal=diffraction_setup_dabax.normalSurface(),
                                                    bragg_angle=diffraction_setup_dabax.angleBragg(energy),
                                                    psi_0=psi_0,
                                                    psi_H=psi_H,
                                                    psi_H_bar=psi_H_bar,
                                                    thickness=diffraction_setup_dabax.thickness(),
                                                    d_spacing=diffraction_setup_dabax.dSpacing() * 1e-10)

        complex_amplitudes = perfect_crystal.calculateDiffraction(photon, method=calculation_method)

        deviations[ia] = deviation
        intensityS_dabax[ia] = complex_amplitudes['S'].intensity()  # 0.0 # coeffs_dabax['S'].intensity()
        intensityP_dabax[ia] = complex_amplitudes['P'].intensity()  # 0.0 # coeffs_dabax['P'].intensity()




    # plot results
    import matplotlib.pylab as plt
    plt.plot(1e6*deviations,intensityS_dabax)
    plt.plot(1e6*deviations,intensityP_dabax)
    plt.xlabel("deviation angle [urad]")
    plt.ylabel("Reflectivity")
    plt.legend(["Sigma-polarization DABAX","Pi-polarization DABAX"])
    plt.show()


#
# main
#
if __name__ == "__main__":

    calculation_method = 1 # 0=Zachariasen, 1=Guigay

    calculate_simple_diffraction_angular_scan_accelerated(calculation_method=calculation_method)

