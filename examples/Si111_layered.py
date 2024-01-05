
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



from crystalpy.diffraction.GeometryType import BraggDiffraction, BraggTransmission, LaueDiffraction, LaueTransmission
from crystalpy.diffraction.DiffractionSetupXraylib import DiffractionSetupXraylib
from crystalpy.diffraction.Diffraction import Diffraction

import scipy.constants as codata

from crystalpy.util.Vector import Vector
from crystalpy.util.Photon import Photon

#
def calculate_simple_diffraction(calculation_method=0, calculation_strategy_flag=1):

    # Create a diffraction setup.

    thickness = 2e-6

    print("\nCreating a diffraction setup...")
    diffraction_setup_r = DiffractionSetupXraylib(geometry_type       = BraggDiffraction(),  # GeometryType object
                                               crystal_name           = "Si",                             # string
                                               thickness              = thickness,                             # meters
                                               miller_h               = 1,                                # int
                                               miller_k               = 1,                                # int
                                               miller_l               = 1,                                # int
                                               asymmetry_angle        = 0,#10.0*numpy.pi/180.,            # radians
                                               azimuthal_angle        = 0.0)                              # radians                            # int

    diffraction_setup_r_half = DiffractionSetupXraylib(geometry_type  = BraggDiffraction(),  # GeometryType object
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
    complex_amplitude_half = numpy.zeros(angle_deviation_points, dtype=complex)
    complex_amplitude      = numpy.zeros(angle_deviation_points, dtype=complex)
    complex_amplitude_bis  = numpy.zeros(angle_deviation_points, dtype=complex)
    complex_amplitude_ter  = numpy.zeros(angle_deviation_points, dtype=complex)



    for ia in range(angle_deviation_points):
        deviation = angle_deviation_min + ia * angle_step
        angle = deviation  + bragg_angle

        # calculate the components of the unitary vector of the incident photon scan
        # Note that diffraction plane is YZ
        yy = numpy.cos(angle)
        zz = - numpy.abs(numpy.sin(angle))
        photon = Photon(energy_in_ev=energy,direction_vector=Vector(0.0,yy,zz))

        # perform the calculation
        coeffs_r_half = Diffraction.calculateDiffractedComplexAmplitudes(diffraction_setup_r_half, photon,
                                                                         calculation_method=calculation_method,
                                                                         is_thick=0,
                                                                         use_transfer_matrix=True,
                                                                         calculation_strategy_flag=calculation_strategy_flag)


        coeffs_r      = Diffraction.calculateDiffractedComplexAmplitudes(diffraction_setup_r, photon,
                                                                         calculation_method=calculation_method,
                                                                         is_thick=0,
                                                                         use_transfer_matrix=True,
                                                                         calculation_strategy_flag=calculation_strategy_flag)

        # 0    1     2    3
        # 11   12    21   22
        # t    t_bar r    r_bar

        complex_amplitude_half[ia] = coeffs_r_half['s21_s']
        complex_amplitude[ia]      = coeffs_r['s21_s']

        #
        # calculate the new transfer matrix and the scattering matrix from it
        #

        # retrieves transfer matrix
        m11 = coeffs_r_half['m11_s']
        m12 = coeffs_r_half['m12_s']
        m21 = coeffs_r_half['m21_s']
        m22 = coeffs_r_half['m22_s']

        # squares it
        m11_bis = (m11**2 + m12 * m21)
        m12_bis = (m11 * m12 + m12 * m22)
        m21_bis = (m21 * m11 + m22 * m21)
        m22_bis = (m21 * m12 + m22**2)

        # just a check
        if True:
            assert (numpy.abs(coeffs_r['m11_s'] - m11_bis) < 1e-10)
            assert (numpy.abs(coeffs_r['m12_s'] - m12_bis) < 1e-10)
            assert (numpy.abs(coeffs_r['m21_s'] - m21_bis) < 1e-10)
            assert (numpy.abs(coeffs_r['m22_s'] - m22_bis) < 1e-10)

        # calculate scattering matrix from transfer matrix (eq 30 in Guigay and Sanchez del Rio)
        from crystalpy.diffraction.PerfectCrystalDiffraction import PerfectCrystalDiffraction
        S11, S12, S21, S22 = PerfectCrystalDiffraction.calculateScatteringMatrixFromTransferMatrix((m11_bis, m12_bis, m21_bis, m22_bis))

        # store complex amplitude of reflectivity
        complex_amplitude_bis[ia] = S21



        #
        # using infinite series from scattering matrix
        #
        t = coeffs_r_half['s11_s']
        t_bar = coeffs_r_half['s22_s']
        r = coeffs_r_half['s21_s']
        r_bar = coeffs_r_half['s12_s']

        complex_amplitude_ter[ia] = r * (1 + t * t_bar / (1 - r * r_bar))

        deviations[ia] = deviation

    # plot results
    from srxraylib.plot.gol import plot

    # plot(1e6 * deviations, numpy.abs(complex_amplitude_half) ** 2,
    #      1e6 * deviations, numpy.abs(complex_amplitude) ** 2,
    #      1e6 * deviations, numpy.abs(complex_amplitude_bis) ** 2,
    #      1e6 * deviations, numpy.abs(complex_amplitude_ter) ** 2,
    #      legend=['half','single','single by transfer matrix multiplication', 'single by series']
    #     )

    plot(1e6 * deviations, numpy.abs(complex_amplitude_half) ** 2,
         1e6 * deviations, numpy.abs(complex_amplitude) ** 2,
         1e6 * deviations, numpy.abs(complex_amplitude_ter) ** 2,
         legend=['T','2T','2T by series'],
         linestyle=[None,None,''],
         marker=[None,None,'+'],
         xtitle=r'$\theta$-$\theta_B$ [$\mu$ rad]',
         ytitle=r'$|r_2|^2$',
        )


#
# main
#
if __name__ == "__main__":

    calculation_method = 1 # 0=Zachariasen, 1=Guigay
    calculation_strategy_flag = 0 # 0=mpmath 1=numpy 2=numpy-truncated
    calculate_simple_diffraction(calculation_method=calculation_method, calculation_strategy_flag=calculation_strategy_flag)

