
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



from crystalpy.diffraction.GeometryType import LaueDiffraction, LaueTransmission, BraggDiffraction, BraggTransmission
from crystalpy.diffraction.DiffractionSetupXraylib import DiffractionSetupXraylib
from crystalpy.diffraction.Diffraction import Diffraction


from crystalpy.util.Vector import Vector


from crystalpy.util.Photon import Photon




#
def calculate_simple_diffraction(geometry_type=LaueDiffraction(), asymmetry_angle=numpy.radians(65), thickness=10e-6,
                                 calculation_method=1, plot_phase=0, calculation_strategy_flag=0, is_thick=0,
                                 save=''):

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
    angle_deviation_min    = -250e-6                          # radians
    angle_deviation_max    =  250e-6                           # radians
    angle_deviation_points = 500

    angle_step = (angle_deviation_max-angle_deviation_min)/angle_deviation_points

    #
    # gets Bragg angle needed to create deviation's scan
    #
    bragg_angle = diffraction_setup.angleBragg(energy)

    print("Bragg angle for E=%f eV is %f deg"%(energy,bragg_angle*180.0/numpy.pi))


    # Create a Diffraction object (the calculator)
    diffraction = Diffraction()

    # initialize arrays for storing outputs
    deviations = numpy.zeros(angle_deviation_points)
    intensityS = numpy.zeros(angle_deviation_points)
    intensityP = numpy.zeros(angle_deviation_points)
    phaseS = numpy.zeros(angle_deviation_points)
    phaseP = numpy.zeros(angle_deviation_points)

    for ia in range(angle_deviation_points):
        deviation = angle_deviation_min + ia * angle_step
        angle = deviation  + bragg_angle + asymmetry_angle
        # print("\n>>>")
        # print("deviation: ", deviation)
        # print("Alpha approx: ", -2 * deviation * numpy.sin(2 * bragg_angle))

        # calculate the components of the unitary vector of the incident photon scan
        # Note that diffraction plane is YZ
        yy = numpy.cos(angle)
        zz = - numpy.abs(numpy.sin(angle))
        photon = Photon(energy_in_ev=energy,direction_vector=Vector(0.0,yy,zz))

        # perform the calculation
        coeffs = diffraction.calculateDiffractedComplexAmplitudes(diffraction_setup, photon,
                                                                  is_thick=is_thick,
                                                                  calculation_method=calculation_method,
                                                                  calculation_strategy_flag=calculation_strategy_flag)
        # print(">>>> ", coeffs["transfer_matrix_s"])
        # store results
        deviations[ia] = deviation
        intensityS[ia] = numpy.abs(coeffs['S']) ** 2
        intensityP[ia] = numpy.abs(coeffs['P']) ** 2
        phaseS[ia] = numpy.angle( numpy.array(coeffs['S'], dtype=complex) )
        phaseP[ia] = numpy.angle( numpy.array(coeffs['P'], dtype=complex) )


    # plot results
    import matplotlib.pylab as plt
    if geometry_type == LaueDiffraction():
        title = "Reflectance"
    elif geometry_type == LaueTransmission():
            title = "Transmittance"
    elif geometry_type == BraggDiffraction():
        title = "Reflectance"
    elif geometry_type == BraggTransmission():
            title = "Transmittance"

    plt.plot(1e6*deviations,intensityS)
    plt.plot(1e6*deviations,intensityP)
    plt.xlabel(r'$\theta$-$\theta_B$ [$\mu$rad]')
    plt.ylabel(title)
    # plt.title(r'Si 111; E = %g eV; $\alpha$=%g deg; $t_c$=%g $\mu$m' % (energy,numpy.rad2deg(asymmetry_angle),thickness*1e6) )
    plt.legend(["Sigma-polarization","Pi-polarization"])
    if save != "": plt.savefig(save)
    plt.show()

    if plot_phase:
        plt.plot(1e6*deviations,phaseS)
        plt.plot(1e6*deviations,phaseP)
        plt.plot(1e6 * deviations, phaseS - phaseP)
        plt.xlabel(r'$\theta$-$\theta_B$ [$\mu$rad]')
        plt.ylabel("Phase "+title)
        plt.title(r'Si 111; E = %g eV; $\alpha$=%g deg; $t_c$=%g $\mu$m' % (energy,numpy.rad2deg(asymmetry_angle),thickness*1e6) )
        plt.legend(["Sigma-polarization","Pi-polarization", "Sigma-Pi"])
        plt.show()


#
# main
#
if __name__ == "__main__":
    import matplotlib.pylab as plt
    plt.rcParams.update({'font.size': 16})
    plt.rcParams.update({'figure.figsize': (10, 6)})

    calculation_method = 1  # 0=Zachariasen, 1=Guigay
    calculation_strategy_flag = 0  # 0=mpmath 1=numpy 2=numpy-truncated

    calculate_simple_diffraction(geometry_type=LaueDiffraction(),   asymmetry_angle=numpy.radians(65), thickness=10e-6, calculation_method=calculation_method, plot_phase=1, calculation_strategy_flag=calculation_strategy_flag, save="Laue_1.pdf")
    calculate_simple_diffraction(geometry_type=LaueTransmission(),  asymmetry_angle=numpy.radians(65), thickness=10e-6, calculation_method=calculation_method, plot_phase=1, calculation_strategy_flag=calculation_strategy_flag, save="Laue_2.pdf")
    calculate_simple_diffraction(geometry_type=BraggDiffraction(),  asymmetry_angle=numpy.radians( 0), thickness=10e-6, calculation_method=calculation_method, plot_phase=1, calculation_strategy_flag=calculation_strategy_flag, save="Bragg_1.pdf")
    calculate_simple_diffraction(geometry_type=BraggTransmission(), asymmetry_angle=numpy.radians( 0), thickness=10e-6, calculation_method=calculation_method, plot_phase=1, calculation_strategy_flag=calculation_strategy_flag, save="Bragg_2.pdf")

    calculate_simple_diffraction(geometry_type=BraggDiffraction(),  asymmetry_angle=numpy.radians( 0), thickness=10e-6, calculation_method=calculation_method, plot_phase=1, calculation_strategy_flag=calculation_strategy_flag, is_thick=1)