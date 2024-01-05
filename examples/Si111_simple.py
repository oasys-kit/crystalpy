
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
from crystalpy.diffraction.DiffractionSetupXraylib import DiffractionSetupXraylib
from crystalpy.diffraction.Diffraction import Diffraction


from crystalpy.util.Vector import Vector
from crystalpy.util.Photon import Photon




#
def calculate_simple_diffraction(calculation_method=0, calculation_strategy_flag=0):

    # Create a diffraction setup.

    print("\nCreating a diffraction setup...")
    diffraction_setup = DiffractionSetupXraylib(geometry_type          = BraggDiffraction(),  # GeometryType object
                                               crystal_name           = "Si",                             # string
                                               thickness              = 1e-2,                             # meters
                                               miller_h               = 1,                                # int
                                               miller_k               = 1,                                # int
                                               miller_l               = 1,                                # int
                                               asymmetry_angle        = 0,#10.0*numpy.pi/180.,                              # radians
                                               azimuthal_angle        = 0.0)                              # radians                            # int


    energy                 = 8000.0                           # eV
    angle_deviation_min    = -100e-6                          # radians
    angle_deviation_max    = 100e-6                           # radians
    angle_deviation_points = 500

    angle_step = (angle_deviation_max-angle_deviation_min)/angle_deviation_points

    #
    # gets Bragg angle needed to create deviation's scan
    #
    bragg_angle = diffraction_setup.angleBragg(energy)

    print("Bragg angle for E=%f eV is %f deg"%(energy,bragg_angle*180.0/numpy.pi))



    # initialize arrays for storing outputs
    deviations = numpy.zeros(angle_deviation_points)
    intensityS = numpy.zeros(angle_deviation_points)
    intensityP = numpy.zeros(angle_deviation_points)

    for ia in range(angle_deviation_points):
        deviation = angle_deviation_min + ia * angle_step
        angle = deviation  + bragg_angle

        # calculate the components of the unitary vector of the incident photon scan
        # Note that diffraction plane is YZ
        yy = numpy.cos(angle)
        zz = - numpy.abs(numpy.sin(angle))
        photon = Photon(energy_in_ev=energy,direction_vector=Vector(0.0,yy,zz))

        # perform the calculation
        coeffs = Diffraction.calculateDiffractedComplexAmplitudes(diffraction_setup, photon,
                                                                  calculation_method=calculation_method,
                                                                  calculation_strategy_flag=calculation_strategy_flag)

        # store results
        deviations[ia] = deviation
        intensityS[ia] = numpy.abs(coeffs['S']) ** 2
        intensityP[ia] = numpy.abs(coeffs['P']) ** 2

    # plot results
    import matplotlib.pylab as plt
    plt.plot(1e6*deviations,intensityS)
    plt.plot(1e6*deviations,intensityP)
    plt.xlabel("deviation angle [urad]")
    plt.ylabel("Reflectivity")
    plt.legend(["Sigma-polarization","Pi-polarization"])
    plt.show()


#
# main
#
if __name__ == "__main__":

    calculation_method = 1 # 0=Zachariasen, 1=Guigay
    calculation_strategy_flag = 2  # 0=mpmath 1=numpy 2=numpy-truncated
    # calculate_simple_diffraction(calculation_method=calculation_method, calculation_strategy_flag=calculation_strategy_flag)



    import numpy
    from crystalpy.util.calc_xcrystal import calc_xcrystal_angular_scan, calc_xcrystal_energy_scan, calc_xcrystal_alphazachariasen_scan

    bunch_out_dict, diffraction_setup, deviations = calc_xcrystal_angular_scan(
        # material_constants_library_flag=self.material_constants_library_flag,
        crystal_name              = 'Si',
        thickness                 = 0.008, #1e-2, #0.007,
        miller_h                  = 1,
        miller_k                  = 1,
        miller_l                  = 1,
        asymmetry_angle           = 0.0,
        energy                    = 8000.0,
        angle_deviation_min       = -0.0001,
        angle_deviation_max       = 0.0001,
        angle_deviation_points    = 200,
        angle_center_flag         = 2,
        calculation_method        = 1, # 0=Zachariasen, 1=Guigay
        is_thick                  = 0,
        use_transfer_matrix       = 0,
        geometry_type_index       = 0,
        calculation_strategy_flag = 2, # 0=mpmath 1=numpy 2=numpy-truncated
                )

    tmp = numpy.zeros((bunch_out_dict["energies"].size,7))
    tmp[:, 0] = deviations / 1e-06
    tmp[:, 1] = 8000.0
    tmp[:, 2] = bunch_out_dict["phaseP"]
    tmp[:, 3] = bunch_out_dict["phaseS"]
    # tmp[:, 4] = circular polarization
    tmp[:, 5] = bunch_out_dict["intensityP"]
    tmp[:, 6] = bunch_out_dict["intensityS"]

    from srxraylib.plot.gol import plot
    plot(tmp[:,0], tmp[:,6], tmp[:,0], tmp[:,5], xtitle="angle", legend=["S-pol","P-pol"])


