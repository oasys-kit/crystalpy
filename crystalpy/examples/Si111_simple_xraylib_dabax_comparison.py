
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
from crystalpy.diffraction.Diffraction import Diffraction


from crystalpy.util.Vector import Vector
from crystalpy.util.Photon import Photon
from crystalpy.diffraction.PerfectCrystalDiffraction import PerfectCrystalDiffraction


from dabax.dabax_xraylib import DabaxXraylib

import time

#
def calculate_simple_diffraction_angular_scan():

    # Create a diffraction setup.

    print("\nCreating a diffraction setup...")
    diffraction_setup = DiffractionSetup(geometry_type          = BraggDiffraction(),  # GeometryType object
                                               crystal_name           = "Si",                             # string
                                               thickness              = 1e-2,                             # meters
                                               miller_h               = 1,                                # int
                                               miller_k               = 1,                                # int
                                               miller_l               = 1,                                # int
                                               asymmetry_angle        = 0,#10.0*numpy.pi/180.,                              # radians
                                               azimuthal_angle        = 0.0)                              # radians                            # int

    import socket
    if socket.getfqdn().find("esrf") >= 0:
        dx = DabaxXraylib(dabax_repository="http://ftp.esrf.fr/pub/scisoft/DabaxFiles/")
    else:
        dx = DabaxXraylib()

    diffraction_setup_dabax = DiffractionSetupDabax(geometry_type          = BraggDiffraction(),  # GeometryType object
                                               crystal_name           = "Si",                             # string
                                               thickness              = 1e-2,                             # meters
                                               miller_h               = 1,                                # int
                                               miller_k               = 1,                                # int
                                               miller_l               = 1,                                # int
                                               asymmetry_angle        = 0,#10.0*numpy.pi/180.,                              # radians
                                               azimuthal_angle        = 0.0,
                                               dabax=dx)                              # radians


    energy                 = 8000.0                           # eV
    angle_deviation_min    = -100e-6                          # radians
    angle_deviation_max    = 100e-6                           # radians
    angle_deviation_points = 50

    angle_step = (angle_deviation_max-angle_deviation_min)/angle_deviation_points

    #
    # gets Bragg angle needed to create deviation's scan
    #
    bragg_angle = diffraction_setup.angleBragg(energy)

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
        coeffs = diffraction.calculateDiffractedComplexAmplitudes(diffraction_setup, photon)
        coeffs_dabax = diffraction_dabax.calculateDiffractedComplexAmplitudes(diffraction_setup_dabax, photon)

        # store results
        deviations[ia] = deviation
        intensityS[ia] = coeffs['S'].intensity()
        intensityP[ia] = coeffs['P'].intensity()
        intensityS_dabax[ia] = coeffs_dabax['S'].intensity()
        intensityP_dabax[ia] = coeffs_dabax['P'].intensity()

    # plot results
    import matplotlib.pylab as plt
    plt.plot(1e6*deviations,intensityS)
    plt.plot(1e6*deviations,intensityP)
    plt.plot(1e6*deviations,intensityS_dabax)
    plt.plot(1e6*deviations,intensityP_dabax)
    plt.xlabel("deviation angle [urad]")
    plt.ylabel("Reflectivity")
    plt.legend(["Sigma-polarization XRAYLIB","Pi-polarization XRAYLIB","Sigma-polarization DABAX","Pi-polarization DABAX"])
    plt.show()


def calculate_simple_diffraction_angular_scan_accelerated():

    # Create a diffraction setup.

    print("\nCreating a diffraction setup...")
    diffraction_setup = DiffractionSetup(geometry_type          = BraggDiffraction(),  # GeometryType object
                                               crystal_name           = "Si",                             # string
                                               thickness              = 1e-2,                             # meters
                                               miller_h               = 1,                                # int
                                               miller_k               = 1,                                # int
                                               miller_l               = 1,                                # int
                                               asymmetry_angle        = 0,#10.0*numpy.pi/180.,                              # radians
                                               azimuthal_angle        = 0.0)                              # radians                            # int

    import socket
    if socket.getfqdn().find("esrf") >= 0:
        dx = DabaxXraylib(dabax_repository="http://ftp.esrf.fr/pub/scisoft/DabaxFiles/")
    else:
        dx = DabaxXraylib()

    diffraction_setup_dabax = DiffractionSetupDabax(geometry_type          = BraggDiffraction(),  # GeometryType object
                                               crystal_name           = "Si",                             # string
                                               thickness              = 1e-2,                             # meters
                                               miller_h               = 1,                                # int
                                               miller_k               = 1,                                # int
                                               miller_l               = 1,                                # int
                                               asymmetry_angle        = 0,#10.0*numpy.pi/180.,                              # radians
                                               azimuthal_angle        = 0.0,
                                               dabax=dx)                              # radians


    energy                 = 8000.0                           # eV
    angle_deviation_min    = -100e-6                          # radians
    angle_deviation_max    = 100e-6                           # radians
    angle_deviation_points = 50

    angle_step = (angle_deviation_max-angle_deviation_min)/angle_deviation_points

    #
    # gets Bragg angle needed to create deviation's scan
    #
    bragg_angle = diffraction_setup.angleBragg(energy)

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
        coeffs = diffraction.calculateDiffractedComplexAmplitudes(diffraction_setup, photon)

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

        complex_amplitudes = perfect_crystal.calculateDiffraction(photon)

        deviations[ia] = deviation
        intensityS_dabax[ia] = complex_amplitudes['S'].intensity()  # 0.0 # coeffs_dabax['S'].intensity()
        intensityP_dabax[ia] = complex_amplitudes['P'].intensity()  # 0.0 # coeffs_dabax['P'].intensity()









    # plot results
    import matplotlib.pylab as plt
    plt.plot(1e6*deviations,intensityS)
    plt.plot(1e6*deviations,intensityP)
    plt.plot(1e6*deviations,intensityS_dabax)
    plt.plot(1e6*deviations,intensityP_dabax)
    plt.xlabel("deviation angle [urad]")
    plt.ylabel("Reflectivity")
    plt.legend(["Sigma-polarization XRAYLIB","Pi-polarization XRAYLIB","Sigma-polarization DABAX","Pi-polarization DABAX"])
    plt.show()


def calculate_simple_diffraction_energy_scan():

    # Create a diffraction setup.

    print("\nCreating a diffraction setup...")
    diffraction_setup = DiffractionSetup(geometry_type          = BraggDiffraction(),  # GeometryType object
                                               crystal_name           = "Si",                             # string
                                               thickness              = 1e-2,                             # meters
                                               miller_h               = 1,                                # int
                                               miller_k               = 1,                                # int
                                               miller_l               = 1,                                # int
                                               asymmetry_angle        = 0,#10.0*numpy.pi/180.,                              # radians
                                               azimuthal_angle        = 0.0)                              # radians                            # int

    import socket
    if socket.getfqdn().find("esrf") >= 0:
        dx = DabaxXraylib(dabax_repository="http://ftp.esrf.fr/pub/scisoft/DabaxFiles/")
    else:
        dx = DabaxXraylib()

    diffraction_setup_dabax = DiffractionSetupDabax(geometry_type          = BraggDiffraction(),  # GeometryType object
                                               crystal_name           = "Si",                             # string
                                               thickness              = 1e-2,                             # meters
                                               miller_h               = 1,                                # int
                                               miller_k               = 1,                                # int
                                               miller_l               = 1,                                # int
                                               asymmetry_angle        = 0,#10.0*numpy.pi/180.,                              # radians
                                               azimuthal_angle        = 0.0,
                                               dabax=dx)                              # radians


    energy                 = 8000.0                           # eV


    angle_deviation_min    = -100e-6                          # radians
    angle_deviation_max    = 100e-6                           # radians
    angle_deviation_points = 50

    angle_step = (angle_deviation_max-angle_deviation_min)/angle_deviation_points

    #
    # gets Bragg angle needed to create deviation's scan
    #
    bragg_angle_corrected = diffraction_setup.angleBraggCorrected(energy)

    print("Bragg angle corrected for E=%f eV is %f deg"%(energy,bragg_angle_corrected*180.0/numpy.pi))


    DeltaE = energy * 1e-4

    npoints = 100
    energies = numpy.linspace(energy-3*DeltaE, energy+3*DeltaE, npoints)

    # Create a Diffraction object (the calculator)
    diffraction = Diffraction()
    diffraction_dabax = Diffraction()

    # initialize arrays for storing outputs
    intensityS =       numpy.zeros(npoints)
    intensityP =       numpy.zeros(npoints)
    intensityS_dabax = numpy.zeros(npoints)
    intensityP_dabax = numpy.zeros(npoints)


    t0 = time.time()
    for ia in range(npoints):

        # calculate the components of the unitary vector of the incident photon scan
        # Note that diffraction plane is YZ
        yy = numpy.cos(bragg_angle_corrected)
        zz = - numpy.abs(numpy.sin(bragg_angle_corrected))
        photon = Photon(energy_in_ev=energies[ia],direction_vector=Vector(0.0,yy,zz))

        # perform the calculation
        coeffs = diffraction.calculateDiffractedComplexAmplitudes(diffraction_setup, photon)

        # store results
        intensityS[ia] = coeffs['S'].intensity()
        intensityP[ia] = coeffs['P'].intensity()


    t1 = time.time()
    for ia in range(npoints):

        # calculate the components of the unitary vector of the incident photon scan
        # Note that diffraction plane is YZ
        yy = numpy.cos(bragg_angle_corrected)
        zz = - numpy.abs(numpy.sin(bragg_angle_corrected))
        photon = Photon(energy_in_ev=energies[ia],direction_vector=Vector(0.0,yy,zz))

        # perform the calculation
        coeffs_dabax = diffraction_dabax.calculateDiffractedComplexAmplitudes(diffraction_setup_dabax, photon)

        # store results
        intensityS_dabax[ia] = coeffs_dabax['S'].intensity()
        intensityP_dabax[ia] = coeffs_dabax['P'].intensity()

    t2 = time.time()

    # plot results
    import matplotlib.pylab as plt
    plt.plot(energies,intensityS)
    plt.plot(energies,intensityP)
    plt.plot(energies,intensityS_dabax)
    plt.plot(energies,intensityP_dabax)
    plt.xlabel("photon energy [eV]")
    plt.ylabel("Reflectivity")
    plt.legend(["Sigma-polarization XRAYLIB","Pi-polarization XRAYLIB","Sigma-polarization DABAX","Pi-polarization DABAX"])
    plt.show()

    print("Total time, Time per points XRAYLIB: ", t1-t0, (t1-t0) / npoints)
    print("Total time, Time per points DABAX: ", t2-t1, (t2-t1) / npoints)

def calculate_simple_diffraction_energy_scan_accelerated():

    # Create a diffraction setup.

    print("\nCreating a diffraction setup...")
    diffraction_setup = DiffractionSetup(geometry_type          = BraggDiffraction(),  # GeometryType object
                                               crystal_name           = "Si",                             # string
                                               thickness              = 1e-2,                             # meters
                                               miller_h               = 1,                                # int
                                               miller_k               = 1,                                # int
                                               miller_l               = 1,                                # int
                                               asymmetry_angle        = 0,#10.0*numpy.pi/180.,                              # radians
                                               azimuthal_angle        = 0.0)                              # radians                            # int

    import socket
    if socket.getfqdn().find("esrf") >= 0:
        dx = DabaxXraylib(dabax_repository="http://ftp.esrf.fr/pub/scisoft/DabaxFiles/")
    else:
        dx = DabaxXraylib()

    diffraction_setup_dabax = DiffractionSetupDabax(geometry_type          = BraggDiffraction(),  # GeometryType object
                                               crystal_name           = "Si",                             # string
                                               thickness              = 1e-2,                             # meters
                                               miller_h               = 1,                                # int
                                               miller_k               = 1,                                # int
                                               miller_l               = 1,                                # int
                                               asymmetry_angle        = 0,#10.0*numpy.pi/180.,                              # radians
                                               azimuthal_angle        = 0.0,
                                               dabax=dx)                              # radians


    energy                 = 8000.0                           # eV


    angle_deviation_min    = -100e-6                          # radians
    angle_deviation_max    = 100e-6                           # radians
    angle_deviation_points = 50

    angle_step = (angle_deviation_max-angle_deviation_min)/angle_deviation_points

    #
    # gets Bragg angle needed to create deviation's scan
    #
    bragg_angle_corrected = diffraction_setup.angleBraggCorrected(energy)

    print("Bragg angle corrected for E=%f eV is %f deg"%(energy,bragg_angle_corrected*180.0/numpy.pi))


    DeltaE = energy * 1e-4

    npoints = 100
    energies = numpy.linspace(energy-3*DeltaE, energy+3*DeltaE, npoints)

    # Create a Diffraction object (the calculator)
    diffraction = Diffraction()
    diffraction_dabax = Diffraction()

    # initialize arrays for storing outputs
    intensityS =       numpy.zeros(npoints)
    intensityP =       numpy.zeros(npoints)
    intensityS_dabax = numpy.zeros(npoints)
    intensityP_dabax = numpy.zeros(npoints)


    t0 = time.time()
    for ia in range(npoints):

        # calculate the components of the unitary vector of the incident photon scan
        # Note that diffraction plane is YZ
        yy = numpy.cos(bragg_angle_corrected)
        zz = - numpy.abs(numpy.sin(bragg_angle_corrected))
        photon = Photon(energy_in_ev=energies[ia],direction_vector=Vector(0.0,yy,zz))

        # perform the calculation
        coeffs = diffraction.calculateDiffractedComplexAmplitudes(diffraction_setup, photon)

        # store results
        intensityS[ia] = coeffs['S'].intensity()
        intensityP[ia] = coeffs['P'].intensity()



    COOR = []
    ENER = []
    for ia in range(npoints):

        # calculate the components of the unitary vector of the incident photon scan
        # Note that diffraction plane is YZ
        yy = numpy.cos(bragg_angle_corrected)
        zz = - numpy.abs(numpy.sin(bragg_angle_corrected))

        COOR.append([0.0,yy,zz])
        ENER.append(energies[ia])


    t1 = time.time()

    Psi_0, Psi_H, Psi_H_bar = diffraction_setup_dabax.psiAll(ENER)

    for ia in range(npoints):
        # perform the calculation
        incoming_photon = Photon(energy_in_ev=ENER[ia],direction_vector=Vector(COOR[ia][0],COOR[ia][1],COOR[ia][2]))
        energy = ENER[ia]

        # psi_0, psi_H, psi_H_bar = diffraction_setup_dabax.psiAll(energy)
        psi_0, psi_H, psi_H_bar = Psi_0[ia], Psi_H[ia], Psi_H_bar[ia]

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

        complex_amplitudes = perfect_crystal.calculateDiffraction(incoming_photon)

        intensityS_dabax[ia] = complex_amplitudes['S'].intensity() # 0.0 # coeffs_dabax['S'].intensity()
        intensityP_dabax[ia] = complex_amplitudes['P'].intensity() # 0.0 # coeffs_dabax['P'].intensity()

    t2 = time.time()

    # plot results
    import matplotlib.pylab as plt
    plt.plot(energies,intensityS)
    plt.plot(energies,intensityP)
    plt.plot(energies,intensityS_dabax)
    plt.plot(energies,intensityP_dabax)
    plt.xlabel("photon energy [eV]")
    plt.ylabel("Reflectivity")
    plt.legend(["Sigma-polarization XRAYLIB","Pi-polarization XRAYLIB","Sigma-polarization DABAX","Pi-polarization DABAX"])
    plt.show()

    print("Total time, Time per points XRAYLIB: ", t1-t0, (t1-t0) / npoints)
    print("Total time, Time per points DABAX: ", t2-t1, (t2-t1) / npoints)


#
# main
#
if __name__ == "__main__":


    # calculate_simple_diffraction_angular_scan()
    # calculate_simple_diffraction_energy_scan()

    calculate_simple_diffraction_angular_scan_accelerated()
    calculate_simple_diffraction_energy_scan_accelerated()

