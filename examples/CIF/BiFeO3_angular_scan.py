
#
#
# This example shows the 204 diffraction by a BiFeO3 crystal at 17 keV:
#
import numpy
import time

from DiffractionSetupCIF import DiffractionSetupCIFdabax, DiffractionSetupCIFxraylib

from crystalpy.diffraction.GeometryType import BraggDiffraction
from crystalpy.diffraction.Diffraction import Diffraction
from crystalpy.util.Vector import Vector
from crystalpy.util.Photon import Photon
from crystalpy.util.ComplexAmplitudePhotonBunch import ComplexAmplitudePhotonBunch
from crystalpy.util.ComplexAmplitudePhoton import ComplexAmplitudePhoton

#
def calculate_simple_diffraction(
        calculation_method=0,  # 0=Zachariasen, 1=Guigay
        ):

    # Create a diffraction setup.

    print("\nCreating a diffraction setup...")
    diffraction_setup = DiffractionSetupCIFdabax(geometry_type          = BraggDiffraction(),  # GeometryType object
                                               crystal_name           = "BiFeO3.cif",                     # string
                                               thickness              = 1e-2,                             # meters
                                               miller_h               = 2,                                # int
                                               miller_k               = 0,                                # int
                                               miller_l               = -4,                               # int
                                               asymmetry_angle        = 0,                                # radians
                                               azimuthal_angle        = 0.0)                              # radians                            # int


    energy                 = 17000.0    # eV
    angle_deviation_min    = -200e-6    # radians
    angle_deviation_max    = 200e-6     # radians
    angle_deviation_points = 5000

    angle_step = (angle_deviation_max-angle_deviation_min) / angle_deviation_points

    #
    # gets Bragg angle needed to create deviation's scan
    #
    bragg_angle = diffraction_setup.angleBragg(energy)

    print("Bragg angle for E=%f eV is %f deg" % (energy, bragg_angle * 180.0 / numpy.pi))

    deviations = numpy.zeros(angle_deviation_points)


    t0 = time.time()
    bunch_in = ComplexAmplitudePhotonBunch()
    K0 = diffraction_setup.vectorK0(energy)
    K0unitary = K0.getNormalizedVector()
    for ia in range(angle_deviation_points):

        deviation = angle_deviation_min + ia * angle_step
        # minus sign in angle is to perform cw rotation when deviation increses
        Vin = K0unitary.rotateAroundAxis(Vector(1, 0, 0), -deviation)
        photon = ComplexAmplitudePhoton(energy_in_ev=energy, direction_vector=Vin)

        bunch_in.addPhoton(photon)
        deviations[ia] = angle_deviation_min + ia * angle_step
    ti1 = time.time()

    coeffs = Diffraction.calculateDiffractedComplexAmplitudes(diffraction_setup,
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

    return deviations, intensityS, intensityP




def calculate_fwhm(x, y):
    y = y / numpy.max(y)  # normalize
    half_max = 0.5
    indices = numpy.where(y >= half_max)[0]

    if len(indices) < 2:
        return 0.0  # No proper peak

    fwhm = x[indices[-1]] - x[indices[0]]
    return fwhm

#
# main
#
if __name__ == "__main__":
    deviations, intensityS, intensityP = calculate_simple_diffraction()

    # plot
    title= "FWHM S-reflectivity: %.3f urad" % (calculate_fwhm(1e6 * deviations, intensityS))
    from srxraylib.plot.gol import plot
    plot(1e6 * deviations, intensityS,
         1e6 * deviations, intensityP,
         xtitle="deviation angle [urad]", ytitle="Reflectivity", title=title,
         legend=["Sigma-polarization","Pi-polarization"],
         grid=1,
         )
