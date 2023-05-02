import numpy


def calc_xcrystal_angular_scan(
        material_constants_library_flag=0,
        crystal_name="Si",
        thickness=1e-2,
        miller_h=1,
        miller_k=1,
        miller_l=1,
        asymmetry_angle=0.0,
        energy=8000.0,
        angle_deviation_min=-100e-6,
        angle_deviation_max=100e-6,
        angle_deviation_points=200,
        calculation_method=0,
        use_corrected_bragg_angle=0,
        do_plot=0,
):
    # Create a diffraction setup.
    print("\nCreating a diffraction setup...")

    print("material_constants_library_flag    = ", material_constants_library_flag)
    print("crystal_name                       = ", crystal_name)
    print("thickness                          = ", thickness)
    print("miller_h                           = ", miller_h)
    print("miller_k                           = ", miller_k)
    print("miller_l                           = ", miller_l)
    print("asymmetry_angle                    = ", asymmetry_angle)
    print("energy                             = ", energy)
    print("angle_deviation_min                = ", angle_deviation_min)
    print("angle_deviation_max                = ", angle_deviation_max)
    print("angle_deviation_points             = ", angle_deviation_points)
    print("calculation_method                 = ", calculation_method)
    print("use_corrected_bragg_angle          = ", use_corrected_bragg_angle)



    # Create a diffraction setup.

    print("\nCreating a diffraction setup...")
    # from crystalpy.diffraction.GeometryType import BraggDiffraction
    # if material_constants_library_flag == 0:

    from crystalpy.diffraction.GeometryType import BraggDiffraction
    from crystalpy.diffraction.DiffractionSetup import DiffractionSetup
    from crystalpy.diffraction.DiffractionSetupDabax import DiffractionSetupDabax
    from crystalpy.diffraction.Diffraction1 import Diffraction1 as Diffraction
    from crystalpy.util.ComplexAmplitudePhoton import ComplexAmplitudePhoton
    from crystalpy.util.ComplexAmplitudePhotonBunch import ComplexAmplitudePhotonBunch
    from crystalpy.util.Vector import Vector

    from srxraylib.plot.gol import plot

    if material_constants_library_flag == 0:
        diffraction_setup = DiffractionSetup(geometry_type          = BraggDiffraction(),       # GeometryType object
                                                   crystal_name           = crystal_name,       # string
                                                   thickness              = thickness,          # meters
                                                   miller_h               = miller_h,           # int
                                                   miller_k               = miller_k,           # int
                                                   miller_l               = miller_l,           # int
                                                   asymmetry_angle        = asymmetry_angle,    # radians
                                                   azimuthal_angle        = 0.0)                # radians                            # int
    else:
        diffraction_setup = DiffractionSetupDabax(geometry_type          = BraggDiffraction(),  # GeometryType object
                                                   crystal_name           = crystal_name,       # string
                                                   thickness              = thickness,          # meters
                                                   miller_h               = miller_h,           # int
                                                   miller_k               = miller_k,           # int
                                                   miller_l               = miller_l,           # int
                                                   asymmetry_angle        = asymmetry_angle,    # radians
                                                   azimuthal_angle        = 0.0)                # radians


    angle_step = (angle_deviation_max-angle_deviation_min)/angle_deviation_points

    bragg_angle = diffraction_setup.angleBragg(energy)

    print("Bragg angle for E=%f eV is %f deg"%(energy,bragg_angle*180.0/numpy.pi))


    # bragg_angle_corrected = diffraction_setup.angleBraggCorrected(energy)
    # print("Bragg angle corrected for E=%f eV is %f deg" % (energy, numpy.degrees(bragg_angle_corrected)))
    # delta_corrected = bragg_angle_corrected - bragg_angle
    # print("delta_corrected: ", delta_corrected)



    # Create a Diffraction object.
    diffraction = Diffraction()


    deviations = numpy.arange(angle_deviation_min, angle_deviation_max, angle_step)


    if True:
        bunch_in = ComplexAmplitudePhotonBunch()
        for ia, deviation in enumerate(deviations):
            # #
            # # get wavevector with incident direction matching Bragg angle
            # #
            # K0 = diffraction_setup.getK0(energy)
            # K0unitary = K0.getNormalizedVector()
            # print("K0unitary", K0.components())
            # # minus sign in angle is to perform cw rotation when deviation increses
            # Vin = K0unitary.rotateAroundAxis(Vector(1,0,0),-deviation)
            # photon = ComplexAmplitudePhoton(energy_in_ev=energy,
            #                                 direction_vector=Vin,
            #                                 Esigma=1,
            #                                 Epi=1)


            photon = ComplexAmplitudePhoton(energy_in_ev=energy,
                                            direction_vector=diffraction_setup.incomingPhotonDirection(energy,deviation,
                                                                use_corrected_bragg_angle=use_corrected_bragg_angle),
                                            Esigma=1,
                                            Epi=1)

            bunch_in.addPhoton( photon )

        bunch_out = diffraction.calculateDiffractedComplexAmplitudePhotonBunch(diffraction_setup,
                                                                               bunch_in,
                                                                               method=calculation_method)

    bunch_out_dict = bunch_out.toDictionary()
    print(bunch_out_dict.keys())

    if do_plot:
        plot(1e6*deviations,bunch_out_dict["intensity"],
             1e6*deviations,bunch_out_dict["intensityS"],
             1e6*deviations,bunch_out_dict["intensityP"],
             xtitle="theta - thetaB [urad]",title="Reflectivity calculation using ComplexAmplitudePhoton method:%d"%calculation_method,
             legend=["Total","Sigma","Pi"])

    return bunch_out_dict, diffraction_setup, deviations



if __name__ == "__main__":
    from srxraylib.plot.gol import set_qt
    set_qt()

    calc_xcrystal_angular_scan(material_constants_library_flag=0, do_plot=True)