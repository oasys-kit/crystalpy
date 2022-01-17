import numpy

from crystalpy.diffraction.GeometryType import BraggDiffraction, LaueDiffraction
from crystalpy.diffraction.DiffractionSetup import DiffractionSetup
from crystalpy.diffraction.Diffraction import Diffraction

import scipy.constants as codata

from crystalpy.util.Vector import Vector

from crystalpy.util.Photon import Photon

def calculate_bragg_dcm(energy_setup=8000.0,energies=numpy.linspace(7990,8010,200)):

    diffraction_setup_r = DiffractionSetup(geometry_type=BraggDiffraction(),  # GeometryType object
                                           crystal_name="Si",  # string
                                           thickness=100e-6,  # meters
                                           miller_h=1,  # int
                                           miller_k=1,  # int
                                           miller_l=1,  # int
                                           asymmetry_angle=0,  # 10.0*numpy.pi/180.,            # radians
                                           azimuthal_angle=0.0)  # radians                            # int

    diffraction = Diffraction()

    scan = numpy.zeros_like(energies)
    r = numpy.zeros_like(energies)


    for i in range(energies.size):
        #
        # gets Bragg angle needed to create deviation's scan
        #
        energy = energies[i]
        bragg_angle = diffraction_setup_r.angleBragg(energy_setup)
        print("Bragg angle for E=%f eV is %f deg" % (energy, bragg_angle * 180.0 / numpy.pi))
        # Create a Diffraction object (the calculator)

        deviation = 0.0 # angle_deviation_min + ia * angle_step
        angle = deviation  + bragg_angle

        # calculate the components of the unitary vector of the incident photon scan
        # Note that diffraction plane is YZ
        yy = numpy.cos(angle)
        zz = - numpy.abs(numpy.sin(angle))
        photon = Photon(energy_in_ev=energy,direction_vector=Vector(0.0,yy,zz))

        # perform the calculation
        coeffs_r = diffraction.calculateDiffractedComplexAmplitudes(diffraction_setup_r, photon)

        scan[i] = energy
        r[i] = numpy.abs( coeffs_r['S'].complexAmplitude() )**4

    return scan,r

def calculate_laue_monochromator(energy_setup=8000.0,energies=numpy.linspace(7990,8010,200)):

    diffraction_setup_r = DiffractionSetup(geometry_type=LaueDiffraction(),  # GeometryType object
                                           crystal_name="Si",  # string
                                           thickness=10e-6,  # meters
                                           miller_h=1,  # int
                                           miller_k=1,  # int
                                           miller_l=1,  # int
                                           asymmetry_angle=numpy.pi/2,  # 10.0*numpy.pi/180.,            # radians
                                           azimuthal_angle=0.0)  # radians                            # int

    diffraction = Diffraction()

    scan = numpy.zeros_like(energies)
    r = numpy.zeros_like(energies)


    for i in range(energies.size):
        #
        # gets Bragg angle needed to create deviation's scan
        #
        energy = energies[i]
        bragg_angle = diffraction_setup_r.angleBragg(energy_setup)
        print("Bragg angle for E=%f eV is %f deg" % (energy, bragg_angle * 180.0 / numpy.pi))
        # Create a Diffraction object (the calculator)

        deviation = 0.0  # angle_deviation_min + ia * angle_step
        angle = deviation + numpy.pi / 2 + bragg_angle

        # calculate the components of the unitary vector of the incident photon scan
        # Note that diffraction plane is YZ
        yy = numpy.cos(angle)
        zz = - numpy.abs(numpy.sin(angle))
        photon = Photon(energy_in_ev=energy,direction_vector=Vector(0.0,yy,zz))

        # perform the calculation
        coeffs_r = diffraction.calculateDiffractedComplexAmplitudes(diffraction_setup_r, photon)

        scan[i] = energy
        r[i] = numpy.abs( coeffs_r['S'].complexAmplitude() )**4

    return scan,r

if __name__ == "__main__":

    from srxraylib.plot.gol import plot

    # scan,r = calculate_bragg_dcm()
    # plot(scan,r, title="Bragg")

    scan,r = calculate_laue_monochromator()
    plot(scan,r, title="Laue")