
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
import time


from crystalpy.diffraction.GeometryType import BraggDiffraction
from crystalpy.diffraction.DiffractionSetupXraylib import DiffractionSetupXraylib
from crystalpy.diffraction.Diffraction import Diffraction


from crystalpy.util.Vector import Vector
from crystalpy.util.Photon import Photon


from crystalpy.util.ComplexAmplitudePhotonBunch import ComplexAmplitudePhotonBunch
from crystalpy.util.ComplexAmplitudePhoton import ComplexAmplitudePhoton

#
def calculate_simple_diffraction(method=0, calculation_method=0):

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
    angle_deviation_points = 5000

    angle_step = (angle_deviation_max-angle_deviation_min)/angle_deviation_points

    #
    # gets Bragg angle needed to create deviation's scan
    #
    bragg_angle = diffraction_setup.angleBragg(energy)

    print("Bragg angle for E=%f eV is %f deg"%(energy,bragg_angle*180.0/numpy.pi))


    # Create a Diffraction object (the calculator)
    diffraction = Diffraction()

    deviations = numpy.zeros(angle_deviation_points)


    if method == 0:
        t0 = time.time()
        intensityS = numpy.zeros(angle_deviation_points)
        intensityP = numpy.zeros(angle_deviation_points)
        # initialize arrays for storing outputs

        for ia in range(angle_deviation_points):
            deviation = angle_deviation_min + ia * angle_step
            angle = deviation  + bragg_angle

            # calculate the components of the unitary vector of the incident photon scan
            # Note that diffraction plane is YZ
            yy = numpy.cos(angle)
            zz = - numpy.abs(numpy.sin(angle))
            photon = Photon(energy_in_ev=energy,direction_vector=Vector(0.0,yy,zz))

            # perform the calculation
            coeffs = diffraction.calculateDiffractedComplexAmplitudes(diffraction_setup, photon, calculation_method=calculation_method)

            # store results
            deviations[ia] = deviation
            intensityS[ia] = numpy.abs(coeffs['S']) ** 2
            intensityP[ia] = numpy.abs(coeffs['P']) ** 2
        t1 = time.time()
    else:
        t0 = time.time()
        bunch_in = ComplexAmplitudePhotonBunch()
        K0 = diffraction_setup.vectorK0(energy)
        K0unitary = K0.getNormalizedVector()
        for ia in range(angle_deviation_points):

            deviation = angle_deviation_min + ia * angle_step

            # angle = deviation + bragg_angle
            # yy = numpy.cos(angle)
            # zz = - numpy.abs(numpy.sin(angle))
            # photon = ComplexAmplitudePhoton(energy_in_ev=energy,direction_vector=Vector(0.0,yy,zz))

            # minus sign in angle is to perform cw rotation when deviation increses
            Vin = K0unitary.rotateAroundAxis(Vector(1, 0, 0), -deviation)
            photon = ComplexAmplitudePhoton(energy_in_ev=energy, direction_vector=Vin)

            bunch_in.addPhoton(photon)
            deviations[ia] = angle_deviation_min + ia * angle_step
        ti1 = time.time()

        if method == 1:
            bunch_out = diffraction.calculateDiffractedComplexAmplitudePhotonBunch(diffraction_setup,
                                                                                   bunch_in,
                                                                                   calculation_method=calculation_method)
            bunch_out_dict = bunch_out.toDictionary()
            intensityS = bunch_out_dict["intensityS"]
            intensityP = bunch_out_dict["intensityP"]
        elif method == 2:
            coeffs = diffraction.calculateDiffractedComplexAmplitudes(diffraction_setup,
                                                                                   bunch_in,
                                                                                   calculation_method=calculation_method)
            intensityS = numpy.abs(coeffs["S"])**2
            intensityP = numpy.abs(coeffs["P"])**2
        t1 = time.time()


    print(">>> spent time: ", t1-t0)
    try:
        print(">>> spent time pure calcul: ", t1 - ti1)
    except:
        pass
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

    method = 2


    calculation_method = 0 # 0=Zachariasen, 1=Guigay
    calculate_simple_diffraction(method=method, calculation_method=calculation_method)

    calculation_method = 1 # 0=Zachariasen, 1=Guigay
    calculate_simple_diffraction(method=method, calculation_method=calculation_method)

