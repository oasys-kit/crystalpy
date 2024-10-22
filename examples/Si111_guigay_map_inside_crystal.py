
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
    angle_deviation_min    = -250e-6                          # radians
    angle_deviation_max    =  250e-6                           # radians
    angle_deviation_points = 500

    s_ratio_deviation_points = 50

    angle_step = (angle_deviation_max-angle_deviation_min)/angle_deviation_points

    #
    # gets Bragg angle needed to create deviation's scan
    #
    bragg_angle = diffraction_setup.angleBragg(energy)

    print("Bragg angle for E=%f eV is %f deg"%(energy,bragg_angle*180.0/numpy.pi))


    # initialize arrays for storing outputs
    deviations = numpy.zeros(angle_deviation_points)

    intensityS = numpy.zeros((angle_deviation_points, s_ratio_deviation_points))
    intensityP = numpy.zeros((angle_deviation_points, s_ratio_deviation_points))

    s_ratios = numpy.linspace(-1.0,0.0,s_ratio_deviation_points)

    for isr in range(s_ratio_deviation_points):
        s_ratio = -s_ratios[isr]

        for ia in range(angle_deviation_points):
            deviation = angle_deviation_min + ia * angle_step
            angle = deviation  + bragg_angle + asymmetry_angle
            # calculate the components of the unitary vector of the incident photon scan
            # Note that diffraction plane is YZ
            yy = numpy.cos(angle)
            zz = - numpy.abs(numpy.sin(angle))
            photon = Photon(energy_in_ev=energy,direction_vector=Vector(0.0,yy,zz))

            if False:
                # perform the calculation
                coeffs = Diffraction.calculateDiffractedComplexAmplitudes(diffraction_setup, photon,
                                                                          calculation_method=1,
                                                                          calculation_strategy_flag=calculation_strategy_flag)
            else:
                # Get PerfectCrystal instance for the current photon.
                perfect_crystal = Diffraction._perfectCrystalForPhoton(diffraction_setup, photon,
                                                                       calculation_strategy_flag=calculation_strategy_flag)
                coeffs = perfect_crystal.calculateDiffractionGuigay(photon, debug=0, s_ratio=s_ratio)


            # store results
            deviations[ia] = deviation
            intensityS[ia, isr] = numpy.abs(coeffs['S']) ** 2
            intensityP[ia, isr] = numpy.abs(coeffs['P']) ** 2



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

    import matplotlib
    import matplotlib.pylab as plt
    matplotlib.rcParams.update({'font.size': 26})
    plot_image(intensityS, deviations*1e6, s_ratios, xtitle=r"$\theta-\theta_B$ [$\mu$rad]", ytitle=r"-$s$/$T$=-$t$/$t_c$", title=title, aspect='auto', figsize=(12, 10), show=0)
    plt.savefig(filename)
    plot(deviations*1e-6, intensityS[:,iii], title=iii)

#
# main
#
if __name__ == "__main__":
    calculation_strategy_flag = 0  # 0=mpmath 1=numpy 2=numpy-truncated

    calculate_diffraction_map_inside_crystal(geometry_type=BraggDiffraction(),  asymmetry_angle=numpy.radians(0), thickness=50e-6, calculation_strategy_flag=calculation_strategy_flag)
    calculate_diffraction_map_inside_crystal(geometry_type=BraggTransmission(), asymmetry_angle=numpy.radians(0), thickness=50e-6, calculation_strategy_flag=calculation_strategy_flag)
    #Laue case
    calculate_diffraction_map_inside_crystal(geometry_type=LaueDiffraction(),  asymmetry_angle=numpy.radians(90), thickness=50e-6, calculation_strategy_flag=calculation_strategy_flag)
    calculate_diffraction_map_inside_crystal(geometry_type=LaueTransmission(), asymmetry_angle=numpy.radians(90), thickness=50e-6, calculation_strategy_flag=calculation_strategy_flag)