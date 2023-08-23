
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
from crystalpy.diffraction.DiffractionSetupDabax import DiffractionSetupDabax
from crystalpy.diffraction.DiffractionSetupShadowPreprocessorV1 import DiffractionSetupShadowPreprocessorV1
from crystalpy.diffraction.DiffractionSetupShadowPreprocessorV2 import DiffractionSetupShadowPreprocessorV2
from crystalpy.diffraction.Diffraction import Diffraction


from crystalpy.util.Vector import Vector
from crystalpy.util.ComplexAmplitudePhoton import ComplexAmplitudePhoton
from crystalpy.util.ComplexAmplitudePhotonBunch import ComplexAmplitudePhotonBunch

from crystalpy.diffraction.PerfectCrystalDiffraction import PerfectCrystalDiffraction


from dabax.dabax_xraylib import DabaxXraylib

import time

#

def calculate_simple_diffraction_energy_scan_accelerated(diffraction_setup, calculation_method=0):


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

    # initialize arrays for storing outputs
    intensityS =       numpy.zeros(npoints)
    intensityP =       numpy.zeros(npoints)


    t0 = time.time()
    bunch = ComplexAmplitudePhotonBunch()

    for ia in range(npoints):

        # calculate the components of the unitary vector of the incident photon scan
        # Note that diffraction plane is YZ
        yy = numpy.cos(bragg_angle_corrected)
        zz = - numpy.abs(numpy.sin(bragg_angle_corrected))
        photon = ComplexAmplitudePhoton(energy_in_ev=energies[ia],direction_vector=Vector(0.0,yy,zz), Esigma=1, Epi=1)

        bunch.addPhoton(photon)

        # # perform the calculation
        # coeffs = diffraction.calculateDiffractedComplexAmplitudes(diffraction_setup, photon, calculation_method=calculation_method)
        #
        # # store results
        # intensityS[ia] = numpy.abs(coeffs['S']) ** 2
        # intensityP[ia] = numpy.abs(coeffs['P']) ** 2

    out_bunch = diffraction.calculateDiffractedComplexAmplitudePhotonBunch(diffraction_setup, bunch,
                                                                        calculation_method=calculation_method)


    t1 = time.time()
    print(">>>>> time: ",  t1-t0, type(diffraction_setup),)

    intensityS = out_bunch.toDictionary()["intensityS"]
    return energies, intensityS


#
# main
#
if __name__ == "__main__":

    create_preprocessor_files = False
    calculation_method = 0 # 0=Zachariasen, 1=Guigay

    # Create a diffraction setup.

    # print("\nCreating a diffraction setup...")
    diffraction_setup_xraylib = DiffractionSetupXraylib(geometry_type          = BraggDiffraction(),  # GeometryType object
                                               crystal_name           = "Si",                             # string
                                               thickness              = 1e-2,                             # meters
                                               miller_h               = 1,                                # int
                                               miller_k               = 1,                                # int
                                               miller_l               = 1,                                # int
                                               asymmetry_angle        = 0,#10.0*numpy.pi/180.,                              # radians
                                               azimuthal_angle        = 0.0)                              # radians                            # int


    diffraction_setup_dabax = DiffractionSetupDabax(geometry_type          = BraggDiffraction(),  # GeometryType object
                                               crystal_name           = "Si",                             # string
                                               thickness              = 1e-2,                             # meters
                                               miller_h               = 1,                                # int
                                               miller_k               = 1,                                # int
                                               miller_l               = 1,                                # int
                                               asymmetry_angle        = 0,#10.0*numpy.pi/180.,                              # radians
                                               azimuthal_angle        = 0.0,
                                               dabax=DabaxXraylib())                              # radians



    if create_preprocessor_files:
        from xoppylib.crystals.create_bragg_preprocessor_file_v1 import create_bragg_preprocessor_file_v1
        from xoppylib.crystals.create_bragg_preprocessor_file_v2 import create_bragg_preprocessor_file_v2
        try:
            import xraylib
        except:
            print("xraylib not available")
        preprocessor_file = "bragg_v1_backup.dat"
        create_bragg_preprocessor_file_v1(interactive=False,
                                              DESCRIPTOR="Si", H_MILLER_INDEX=1, K_MILLER_INDEX=1, L_MILLER_INDEX=1,
                                              TEMPERATURE_FACTOR=1.0,
                                              E_MIN=5000.0, E_MAX=15000.0, E_STEP=100.0,
                                              SHADOW_FILE=preprocessor_file,
                                              material_constants_library=xraylib)

        preprocessor_file = "bragg_v1.dat"
        create_bragg_preprocessor_file_v1(interactive=False,
                                              DESCRIPTOR="Si", H_MILLER_INDEX=1, K_MILLER_INDEX=1, L_MILLER_INDEX=1,
                                              TEMPERATURE_FACTOR=1.0,
                                              E_MIN=5000.0, E_MAX=15000.0, E_STEP=100.0,
                                              SHADOW_FILE=preprocessor_file,
                                              material_constants_library=xraylib)

        preprocessor_file = "bragg_v2_backup.dat"
        create_bragg_preprocessor_file_v2(interactive=False,
                                              DESCRIPTOR="Si", H_MILLER_INDEX=1, K_MILLER_INDEX=1, L_MILLER_INDEX=1,
                                              TEMPERATURE_FACTOR=1.0,
                                              E_MIN=5000.0, E_MAX=15000.0, E_STEP=100.0,
                                              SHADOW_FILE=preprocessor_file,
                                              material_constants_library=xraylib)

        preprocessor_file = "bragg_v2.dat"
        create_bragg_preprocessor_file_v2(interactive=False,
                                              DESCRIPTOR="Si", H_MILLER_INDEX=1, K_MILLER_INDEX=1, L_MILLER_INDEX=1,
                                              TEMPERATURE_FACTOR=1.0,
                                              E_MIN=5000.0, E_MAX=15000.0, E_STEP=100.0,
                                              SHADOW_FILE=preprocessor_file,
                                              material_constants_library=xraylib)


        preprocessor_file_v1 = "bragg_v1.dat"
        preprocessor_file_v2 = "bragg_v2.dat"
    else:
        print("Shadow preprocessor files not created. Using backup files")
        preprocessor_file_v1 = "bragg_v1_backup.dat"
        preprocessor_file_v2 = "bragg_v2_backup.dat"

    diffraction_setup_v1 = DiffractionSetupShadowPreprocessorV1(
                 geometry_type=BraggDiffraction(),
                 crystal_name="Si", thickness=1e-2,
                 miller_h=1, miller_k=1, miller_l=1,
                 asymmetry_angle=0.0,
                 azimuthal_angle=0.0,
                 preprocessor_file=preprocessor_file_v1)

    diffraction_setup_v2 = DiffractionSetupShadowPreprocessorV2(
                 geometry_type=BraggDiffraction(),
                 crystal_name="Si", thickness=1e-2,
                 miller_h=1, miller_k=1, miller_l=1,
                 asymmetry_angle=0.0,
                 azimuthal_angle=0.0,
                 preprocessor_file=preprocessor_file_v2)
    #
    ex, ix = calculate_simple_diffraction_energy_scan_accelerated(diffraction_setup_xraylib, calculation_method=calculation_method)
    ed, id = calculate_simple_diffraction_energy_scan_accelerated(diffraction_setup_dabax, calculation_method=calculation_method)
    e1, i1 = calculate_simple_diffraction_energy_scan_accelerated(diffraction_setup_v1, calculation_method=calculation_method)
    e2, i2 = calculate_simple_diffraction_energy_scan_accelerated(diffraction_setup_v2, calculation_method=calculation_method)
    from srxraylib.plot.gol import plot
    plot(ex, ix,
         ed, id,
         e1, i1,
         e2, i2,
         legend=["xraylib", "dabax", "v1", "v2"])

