import numpy

# used in get_diffraction_setup
from crystalpy.diffraction.GeometryType import BraggDiffraction, LaueTransmission, LaueDiffraction, BraggTransmission
from crystalpy.diffraction.DiffractionSetupXraylib import DiffractionSetupXraylib
from crystalpy.diffraction.DiffractionSetupDabax import DiffractionSetupDabax

# use in calc_xcrystal*
from crystalpy.diffraction.Diffraction import Diffraction
from crystalpy.util.ComplexAmplitudePhoton import ComplexAmplitudePhoton
from crystalpy.util.ComplexAmplitudePhotonBunch import ComplexAmplitudePhotonBunch

def get_diffraction_setup(
        material_constants_library_flag=0,
        crystal_name="Si",
        thickness=1e-2,
        miller_h=1,
        miller_k=1,
        miller_l=1,
        geometry_type_index=0,
        asymmetry_angle=0.0,
        ):

    # Create a diffraction setup.
    print("\nCreating a diffraction setup...")

    if geometry_type_index == 0:
        geometry_type = BraggDiffraction()
    elif geometry_type_index == 1:
        geometry_type = LaueDiffraction()
    elif geometry_type_index == 2:
        geometry_type = BraggTransmission()
    elif geometry_type_index == 3:
        geometry_type = LaueTransmission()

    if material_constants_library_flag == 0:
        diffraction_setup = DiffractionSetupXraylib(geometry_type                = geometry_type,      # GeometryType object
                                                   crystal_name           = crystal_name,       # string
                                                   thickness              = thickness,          # meters
                                                   miller_h               = miller_h,           # int
                                                   miller_k               = miller_k,           # int
                                                   miller_l               = miller_l,           # int
                                                   asymmetry_angle        = asymmetry_angle,    # radians
                                                   azimuthal_angle        = 0.0)                # radians                            # int
    else:
        diffraction_setup = DiffractionSetupDabax(geometry_type           = geometry_type,      # GeometryType object
                                                   crystal_name           = crystal_name,       # string
                                                   thickness              = thickness,          # meters
                                                   miller_h               = miller_h,           # int
                                                   miller_k               = miller_k,           # int
                                                   miller_l               = miller_l,           # int
                                                   asymmetry_angle        = asymmetry_angle,    # radians
                                                   azimuthal_angle        = 0.0)                # radians

    return diffraction_setup

def calc_xcrystal_angular_scan(
        crystal_name="Si",
        miller_h=1,
        miller_k=1,
        miller_l=1,
        thickness=1e-2,
        asymmetry_angle=0.0,
        geometry_type_index=0,
        material_constants_library_flag=0,
        energy=8000.0,
        angle_deviation_min=-100e-6,
        angle_deviation_max=100e-6,
        angle_deviation_points=200,
        angle_center_flag=2, # 0=Absolute angle, 1=Theta Bragg Corrected, 2=Theta Bragg
        calculation_method=0,
        do_plot=0,
):

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
    print("angle_center_flag                  = ", angle_center_flag)
    print("calculation_method                 = ", calculation_method)
    print("geometry_type_index                = ", geometry_type_index)

    # Create a diffraction setup.
    diffraction_setup = get_diffraction_setup(
        crystal_name=crystal_name,
        miller_h=miller_h,
        miller_k=miller_k,
        miller_l=miller_l,
        thickness=thickness,
        asymmetry_angle=asymmetry_angle,
        geometry_type_index=geometry_type_index,
        material_constants_library_flag=material_constants_library_flag,
        )


    bragg_angle = diffraction_setup.angleBragg(energy)
    print("Bragg angle for E=%f eV is %f deg"%(energy,bragg_angle*180.0/numpy.pi))

    bragg_angle_corrected = diffraction_setup.angleBraggCorrected(energy)
    print("Bragg angle corrected for E=%f eV is %f deg" % (energy, numpy.degrees(bragg_angle_corrected)))

    # Create a Diffraction object.
    diffraction = Diffraction()

    deviations = numpy.linspace(angle_deviation_min, angle_deviation_max, angle_deviation_points)

    bunch_in = ComplexAmplitudePhotonBunch()
    for ia, deviation in enumerate(deviations):

        photon = ComplexAmplitudePhoton(energy_in_ev=energy,
                                        direction_vector=diffraction_setup.incomingPhotonDirection(
                                            energy,deviation,
                                            angle_center_flag=angle_center_flag,
                                        ),
                                        Esigma=1,
                                        Epi=1)

        bunch_in.addPhoton(photon)

    bunch_out = diffraction.calculateDiffractedComplexAmplitudePhotonBunch(diffraction_setup,
                                                                           bunch_in,
                                                                           method=calculation_method)

    bunch_out_dict = bunch_out.toDictionary()

    if do_plot:
        from srxraylib.plot.gol import plot
        plot(1e6*deviations,bunch_out_dict["intensity"],
             1e6*deviations,bunch_out_dict["intensityS"],
             1e6*deviations,bunch_out_dict["intensityP"],
             xtitle="theta - thetaB [urad]",title="Reflectivity calculation using ComplexAmplitudePhoton method:%d"%calculation_method,
             legend=["Total","Sigma","Pi"])

    return bunch_out_dict, diffraction_setup, deviations


def calc_xcrystal_alphazachariasen_scan(
        crystal_name="Si",
        miller_h=1,
        miller_k=1,
        miller_l=1,
        thickness=1e-2,
        asymmetry_angle=0.0,
        geometry_type_index=0,
        material_constants_library_flag=0,
        energy=8000.0,
        angle_deviation_min=-2,
        angle_deviation_max=2,
        angle_deviation_points=200,
        calculation_method=0,
        do_plot=0,
):

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
    print("geometry_type_index                = ", geometry_type_index)

    # Create a diffraction setup.
    diffraction_setup = get_diffraction_setup(
        crystal_name=crystal_name,
        miller_h=miller_h,
        miller_k=miller_k,
        miller_l=miller_l,
        thickness=thickness,
        asymmetry_angle=asymmetry_angle,
        geometry_type_index=geometry_type_index,
        material_constants_library_flag=material_constants_library_flag,
        )


    bragg_angle = diffraction_setup.angleBragg(energy)
    print("Bragg angle for E=%f eV is %f deg"%(energy,bragg_angle*180.0/numpy.pi))

    bragg_angle_corrected = diffraction_setup.angleBraggCorrected(energy)
    print("Bragg angle corrected for E=%f eV is %f deg" % (energy, numpy.degrees(bragg_angle_corrected)))

    darwin_half_width = diffraction_setup.darwinHalfwidth(energy)
    print("Darwin half width for E=%f eV is %f deg" % (energy, numpy.degrees(darwin_half_width[0])))


    # Create a Diffraction object.
    diffraction = Diffraction()

    deviations = numpy.linspace(angle_deviation_min, angle_deviation_max, angle_deviation_points)

    bunch_in = ComplexAmplitudePhotonBunch()
    for ia, deviation in enumerate(deviations):

        photon = ComplexAmplitudePhoton(energy_in_ev=energy,
                                        direction_vector=diffraction_setup.incomingPhotonDirection(
                                            energy,
                                            deviation * darwin_half_width[0],
                                            angle_center_flag=1,
                                        ),
                                        Esigma=1,
                                        Epi=1)

        bunch_in.addPhoton(photon)

    bunch_out = diffraction.calculateDiffractedComplexAmplitudePhotonBunch(diffraction_setup,
                                                                           bunch_in,
                                                                           method=calculation_method)

    bunch_out_dict = bunch_out.toDictionary()

    if do_plot:
        from srxraylib.plot.gol import plot
        plot(deviations,bunch_out_dict["intensity"],
             deviations,bunch_out_dict["intensityS"],
             deviations,bunch_out_dict["intensityP"],
             xtitle="y_zac",title="Reflectivity calculation using ComplexAmplitudePhoton method:%d"%calculation_method,
             legend=["Total","Sigma","Pi"])

    return bunch_out_dict, diffraction_setup, deviations


def calc_xcrystal_energy_scan(
        crystal_name="Si",
        miller_h=1,
        miller_k=1,
        miller_l=1,
        thickness=1e-2,
        asymmetry_angle=0.0,
        geometry_type_index=0,
        material_constants_library_flag=0,
        theta=None,
        energy_min=7990,
        energy_max=8010,
        energy_points=200,
        calculation_method=0,
        do_plot=0,
):

    print("material_constants_library_flag    = ", material_constants_library_flag)
    print("crystal_name                       = ", crystal_name)
    print("thickness                          = ", thickness)
    print("miller_h                           = ", miller_h)
    print("miller_k                           = ", miller_k)
    print("miller_l                           = ", miller_l)
    print("asymmetry_angle                    = ", asymmetry_angle)
    print("theta                              = ", theta)
    print("energy_min                         = ", energy_min)
    print("energy_max                         = ", energy_max)
    print("energy_points                      = ", energy_points)
    print("calculation_method                 = ", calculation_method)
    print("geometry_type_index                = ", geometry_type_index)

    # Create a diffraction setup.
    diffraction_setup = get_diffraction_setup(
        crystal_name=crystal_name,
        miller_h=miller_h,
        miller_k=miller_k,
        miller_l=miller_l,
        thickness=thickness,
        asymmetry_angle=asymmetry_angle,
        geometry_type_index=geometry_type_index,
        material_constants_library_flag=material_constants_library_flag,
        )

    if theta is None:
        energy_mean = 0.5 * (energy_max + energy_min)
        theta = diffraction_setup.angleBragg(energy_mean)
        print("Using theta as the Bragg angle for mean E=%f eV, which is %f deg"%(energy_mean, numpy.degrees(theta)))

    # Create a Diffraction object.
    diffraction = Diffraction()

    energies = numpy.linspace(energy_min, energy_max, energy_points)

    bunch_in = ComplexAmplitudePhotonBunch()
    for ia, energy in enumerate(energies):
        bragg_angle = diffraction_setup.angleBragg(energy)
        deviation = theta - bragg_angle
        photon = ComplexAmplitudePhoton(energy_in_ev=energy,
                                        direction_vector=diffraction_setup.incomingPhotonDirection(
                                            energy,
                                            deviation,
                                            angle_center_flag=2,
                                        ),
                                        Esigma=1,
                                        Epi=1)

        bunch_in.addPhoton(photon)

    bunch_out = diffraction.calculateDiffractedComplexAmplitudePhotonBunch(diffraction_setup,
                                                                           bunch_in,
                                                                           method=calculation_method)

    bunch_out_dict = bunch_out.toDictionary()

    if do_plot:
        from srxraylib.plot.gol import plot
        plot(energies,bunch_out_dict["intensity"],
             energies,bunch_out_dict["intensityS"],
             energies,bunch_out_dict["intensityP"],
             xtitle="photon energy [eV]",title="Reflectivity calculation using ComplexAmplitudePhoton method:%d"%calculation_method,
             legend=["Total","Sigma","Pi"])

    return bunch_out_dict, diffraction_setup, energies

def calc_xcrystal_double_scan(
        material_constants_library_flag=0,
        crystal_name="Si",
        thickness=1e-2,
        miller_h=1,
        miller_k=1,
        miller_l=1,
        asymmetry_angle=0.0,
        energy_min=7990,
        energy_max=8010,
        energy_points=200,
        angle_deviation_min=-100e-6,
        angle_deviation_max=100e-6,
        angle_deviation_points=200,
        angle_center_flag=2, # 0=Absolute angle, 1=Theta Bragg Corrected, 2=Theta Bragg
        calculation_method=0,
        geometry_type_index=0,
        do_plot=0,
):


    print("material_constants_library_flag    = ", material_constants_library_flag)
    print("crystal_name                       = ", crystal_name)
    print("thickness                          = ", thickness)
    print("miller_h                           = ", miller_h)
    print("miller_k                           = ", miller_k)
    print("miller_l                           = ", miller_l)
    print("asymmetry_angle                    = ", asymmetry_angle)
    print("energy_min                         = ", energy_min)
    print("energy_max                         = ", energy_max)
    print("energy_points                      = ", energy_points)
    print("angle_deviation_min                = ", angle_deviation_min)
    print("angle_deviation_max                = ", angle_deviation_max)
    print("angle_deviation_points             = ", angle_deviation_points)
    print("angle_center_flag                  = ", angle_center_flag)
    print("calculation_method                 = ", calculation_method)
    print("geometry_type_index                = ", geometry_type_index)


    # Create a diffraction setup.
    diffraction_setup = get_diffraction_setup(
        crystal_name=crystal_name,
        miller_h=miller_h,
        miller_k=miller_k,
        miller_l=miller_l,
        thickness=thickness,
        asymmetry_angle=asymmetry_angle,
        geometry_type_index=geometry_type_index,
        material_constants_library_flag=material_constants_library_flag,
        )

    if energy_points == 1:
        energy_mean = energy_min
    else:
        energy_mean = 0.5 * (energy_min + energy_max)


    # theta = diffraction_setup.angleBragg(energy_mean)
    # print("Using angular deviation center as Bragg angle for mean E=%f eV, which is %f deg"%(energy_mean, numpy.degrees(theta)))

    bragg_angle = diffraction_setup.angleBragg(energy_mean)
    print("Bragg angle for E=%f eV is %f deg"%(energy_mean, bragg_angle*180.0/numpy.pi))

    bragg_angle_corrected = diffraction_setup.angleBraggCorrected(energy_mean)
    print("Bragg angle corrected for E=%f eV is %f deg" % (energy_mean, numpy.degrees(bragg_angle_corrected)))

    # Create a Diffraction object.
    diffraction = Diffraction()

    deviations = numpy.linspace(angle_deviation_min, angle_deviation_max, angle_deviation_points)
    energies   = numpy.linspace(energy_min, energy_max, energy_points)

    bunch_in = ComplexAmplitudePhotonBunch()
    print("Creating bunch...")
    for ie, energy in enumerate(energies):
        for ia, deviation in enumerate(deviations):
            photon = ComplexAmplitudePhoton(energy_in_ev=energy,
                                            direction_vector=diffraction_setup.incomingPhotonDirection(
                                                energy_mean,
                                                deviation,
                                                angle_center_flag=angle_center_flag,
                                            ),
                                            Esigma=1,
                                            Epi=1)

            bunch_in.addPhoton(photon)
    print("Done.")

    print("Calculating diffraction...")
    bunch_out = diffraction.calculateDiffractedComplexAmplitudePhotonBunch(diffraction_setup,
                                                                           bunch_in,
                                                                           method=calculation_method)

    print("Done.")

    bunch_out_dict = bunch_out.toDictionary()

    if do_plot:
        from srxraylib.plot.gol import plot

        if angle_deviation_points == 1:
            plot(energies,bunch_out_dict["intensity"],
                 energies,bunch_out_dict["intensityS"],
                 energies,bunch_out_dict["intensityP"],
                 xtitle="photon energy [eV]",title="Reflectivity calculation using ComplexAmplitudePhoton method:%d"%calculation_method,
                 legend=["Total","Sigma","Pi"])

        if energy_points == 1:
            plot(1e6*deviations,bunch_out_dict["intensity"],
                 1e6*deviations,bunch_out_dict["intensityS"],
                 1e6*deviations,bunch_out_dict["intensityP"],
                 xtitle="theta - thetaB [urad]",title="Reflectivity calculation using ComplexAmplitudePhoton method:%d"%calculation_method,
                 legend=["Total","Sigma","Pi"])

        if (angle_deviation_points > 1) and (energy_points > 1):
            from srxraylib.plot.gol import plot_image
            image = bunch_out_dict["intensityS"]
            image.shape = (energy_points, angle_deviation_points)
            print(image.shape)
            plot_image(image, energies, deviations*1e6, xtitle="photon energy [eV]", ytitle="angle deviation [urad]", aspect='auto')
    return bunch_out_dict, diffraction_setup, energies, deviations

if __name__ == "__main__":
    from srxraylib.plot.gol import set_qt
    set_qt()

    if True:
        calc_xcrystal_angular_scan(material_constants_library_flag=0, do_plot=True)

        calc_xcrystal_angular_scan(material_constants_library_flag=0, geometry_type_index=1, thickness=10e-6, asymmetry_angle=numpy.radians(90), do_plot=True)

        calc_xcrystal_energy_scan(material_constants_library_flag=0, do_plot=True)

    calc_xcrystal_alphazachariasen_scan(do_plot=1)

    if True:
        calc_xcrystal_double_scan(        material_constants_library_flag=0,
            crystal_name="Si",
            thickness=1e-2,
            miller_h=1,
            miller_k=1,
            miller_l=1,
            asymmetry_angle=0.0,
            energy_min=8000,
            energy_max=8010,
            energy_points=1,
            angle_deviation_min=-100e-6,
            angle_deviation_max=100e-6,
            angle_deviation_points=200,
            angle_center_flag=2, # 0=Absolute angle, 1=Theta Bragg Corrected, 2=Theta Bragg
            calculation_method=0,
            geometry_type_index=0,
            do_plot=1,)

        calc_xcrystal_double_scan(        material_constants_library_flag=0,
            crystal_name="Si",
            thickness=1e-2,
            miller_h=1,
            miller_k=1,
            miller_l=1,
            asymmetry_angle=0.0,
            energy_min=7990,
            energy_max=8010,
            energy_points=100,
            angle_deviation_min=0,
            angle_deviation_max=100e-6,
            angle_deviation_points=1,
            angle_center_flag=2, # 0=Absolute angle, 1=Theta Bragg Corrected, 2=Theta Bragg
            calculation_method=0,
            geometry_type_index=0,
            do_plot=1,)

    if False:
        calc_xcrystal_double_scan(        material_constants_library_flag=0,
            crystal_name="Si",
            thickness=1e-2,
            miller_h=1,
            miller_k=1,
            miller_l=1,
            asymmetry_angle=0.0,
            energy_min=7990,
            energy_max=8010,
            energy_points=200,
            angle_deviation_min=-100e-6,
            angle_deviation_max=100e-6,
            angle_deviation_points=100,
            angle_center_flag=2, # 0=Absolute angle, 1=Theta Bragg Corrected, 2=Theta Bragg
            calculation_method=0,
            geometry_type_index=0,
            do_plot=1,)