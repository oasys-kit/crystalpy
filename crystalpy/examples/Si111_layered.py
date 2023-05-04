
#
#
# This example checks that the reflectivity of a crystal of 2 um crystal can be calculated
# from the results of 1 um crystal. For that
#     - calculate the transfer matrix of a 1 um crystal
#     - calculate the transfer matrix of the 2 um crystal, that is the square of the previous one
#     - calculate the scattering matrix
#     - the reflectivity r = is the term 21 of the scattering matrix
#
#
import numpy



from crystalpy.diffraction.GeometryType import BraggDiffraction, BraggTransmission
from crystalpy.diffraction.DiffractionSetupXraylib import DiffractionSetupXraylib
from crystalpy.diffraction.Diffraction import Diffraction

import scipy.constants as codata

from crystalpy.util.Vector import Vector


from crystalpy.util.Photon import Photon




#
def calculate_simple_diffraction(calculation_method=0):

    # Create a diffraction setup.

    thickness = 2e-6

    print("\nCreating a diffraction setup...")
    diffraction_setup_r = DiffractionSetupXraylib(geometry_type          = BraggDiffraction(),  # GeometryType object
                                               crystal_name           = "Si",                             # string
                                               thickness              = thickness,                             # meters
                                               miller_h               = 1,                                # int
                                               miller_k               = 1,                                # int
                                               miller_l               = 1,                                # int
                                               asymmetry_angle        = 0,#10.0*numpy.pi/180.,            # radians
                                               azimuthal_angle        = 0.0)                              # radians                            # int

    diffraction_setup_r_half = DiffractionSetupXraylib(geometry_type          = BraggDiffraction(),  # GeometryType object
                                               crystal_name           = "Si",                             # string
                                               thickness              = thickness/2,                             # meters
                                               miller_h               = 1,                                # int
                                               miller_k               = 1,                                # int
                                               miller_l               = 1,                                # int
                                               asymmetry_angle        = 0,#10.0*numpy.pi/180.,            # radians
                                               azimuthal_angle        = 0.0)                              # radians                            # int


    energy                 = 8000.0                           # eV
    angle_deviation_min    = -300e-6                          # radians
    angle_deviation_max    = 300e-6                           # radians
    angle_deviation_points = 500

    wavelength = codata.h * codata.c / codata.e / energy

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
    complex_amplitude_half = numpy.zeros(angle_deviation_points)
    complex_amplitude = numpy.zeros(angle_deviation_points)
    complex_amplitude_bis = numpy.zeros(angle_deviation_points)


    for ia in range(angle_deviation_points):
        deviation = angle_deviation_min + ia * angle_step
        angle = deviation  + bragg_angle

        # calculate the components of the unitary vector of the incident photon scan
        # Note that diffraction plane is YZ
        yy = numpy.cos(angle)
        zz = - numpy.abs(numpy.sin(angle))
        photon = Photon(energy_in_ev=energy,direction_vector=Vector(0.0,yy,zz))

        # perform the calculation
        coeffs_r_half = diffraction.calculateDiffractedComplexAmplitudes(diffraction_setup_r_half, photon, method=calculation_method)
        coeffs_r      = diffraction.calculateDiffractedComplexAmplitudes(diffraction_setup_r, photon, method=calculation_method)

        # 0    1     2    3
        # 11   12    21   22
        # t    t_bar r    r_bar
        complex_amplitude_half[ia] = coeffs_r_half['scattering_matrix_s'][2]
        complex_amplitude[ia]      = coeffs_r['scattering_matrix_s'][2]

        # retrieves transfer matrix
        m11, m12, m21, m22 = coeffs_r_half['transfer_matrix_s']

        # squares it
        m11_bis = (m11**2 + m12 * m21)
        m12_bis = (m11 * m12 + m12 * m22)
        m21_bis = (m21 * m11 + m22 * m21)
        m22_bis = (m21 * m12 + m22**2)

        # just a check
        if False:
            M11, M12, M21, M22 =      coeffs_r['transfer_matrix_s']
            print(M11 - m11_bis, M12 - m12_bis,  M21 - m21_bis,   M22 - m22_bis )

        # calculate scattering matrix from transfer matrix (eq 30 in Guigay and Sanchez del Rio)
        from crystalpy.diffraction.PerfectCrystalDiffraction import PerfectCrystalDiffraction
        S11, S12, S21, S22 = PerfectCrystalDiffraction.calculateScatteringMatrixFromTransferMatrix((m11_bis, m12_bis, m21_bis, m22_bis))

        # store complex amplitude of reflectivity
        complex_amplitude_bis[ia] = S21

        deviations[ia] = deviation

    # plot results
    from srxraylib.plot.gol import plot

    plot(1e6 * deviations, numpy.abs(complex_amplitude_half) ** 2,
         1e6 * deviations, numpy.abs(complex_amplitude) ** 2,
         1e6 * deviations, numpy.abs(complex_amplitude_bis) ** 2,
         legend=['half','single','single by matrix multiplication']
        )


#
# main
#
if __name__ == "__main__":

    calculation_method = 1 # 0=Zachariasen, 1=Guigay
    calculate_simple_diffraction(calculation_method=calculation_method)

