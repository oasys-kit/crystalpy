import numpy
import time
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
        dabax=None, # only used if material_constants_library_flag > 0
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
                                                   azimuthal_angle        = 0.0,                # radians
                                                   dabax                  = dabax,
                                                   )

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
        is_thick=0,
        use_transfer_matrix=0,
        do_plot=0,
        calculation_strategy_flag=0,  # 0=mpmath 1=numpy 2=numpy-truncated
        dabax=None, # only used if material_constants_library_flag > 0
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
    print("is_thick                           = ", is_thick)
    print("use_transfer_matrix                = ", use_transfer_matrix)
    print("geometry_type_index                = ", geometry_type_index)
    if dabax is None:
        print("dabax                              = ", dabax)
    else:
        print("dabax                              = \n", dabax.info())

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
        dabax=dabax,
        )


    bragg_angle = diffraction_setup.angleBragg(energy)
    print("Bragg angle for E=%f eV is %f deg"%(energy,bragg_angle*180.0/numpy.pi))

    bragg_angle_corrected = diffraction_setup.angleBraggCorrected(energy)
    print("Bragg angle corrected for E=%f eV is %f deg" % (energy, numpy.degrees(bragg_angle_corrected)))


    deviations = numpy.linspace(angle_deviation_min, angle_deviation_max, angle_deviation_points)

    bunch_in = ComplexAmplitudePhotonBunch()
    for ia, deviation in enumerate(deviations):

        photon = ComplexAmplitudePhoton(energy_in_ev=energy,
                                        direction_vector=diffraction_setup.vectorIncomingPhotonDirection(
                                            energy,deviation,
                                            angle_center_flag=angle_center_flag,
                                        ),
                                        Esigma=1,
                                        Epi=1)

        bunch_in.addPhoton(photon)

    bunch_out = Diffraction.calculateDiffractedComplexAmplitudePhotonBunch(diffraction_setup,
                                                                           bunch_in,
                                                                           calculation_method=calculation_method,
                                                                           is_thick=is_thick,
                                                                           use_transfer_matrix=use_transfer_matrix,
                                                                           calculation_strategy_flag=calculation_strategy_flag,
                                                                           )

    bunch_out_dict = bunch_out.toDictionary()

    if do_plot:
        from srxraylib.plot.gol import plot
        plot(1e6*deviations,bunch_out_dict["intensity"],
             1e6*deviations,bunch_out_dict["intensityS"],
             1e6*deviations,bunch_out_dict["intensityP"],
             xtitle="theta - thetaB [urad]",title="Reflectivity calculation using ComplexAmplitudePhoton calculation_method:%d"%calculation_method,
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
        is_thick=0,
        use_transfer_matrix=0,
        do_plot=0,
        calculation_strategy_flag=0,  # 0=mpmath 1=numpy 2=numpy-truncated
        dabax=None, # only used if material_constants_library_flag > 0
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
    print("is_thick                           = ", is_thick)
    print("use_transfer_matrix                = ", use_transfer_matrix)
    print("geometry_type_index                = ", geometry_type_index)
    if dabax is None:
        print("dabax                              = ", dabax)
    else:
        print("dabax                              = \n", dabax.info())

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
        dabax=dabax,
        )


    bragg_angle = diffraction_setup.angleBragg(energy)
    print("Bragg angle for E=%f eV is %f deg"%(energy,bragg_angle*180.0/numpy.pi))

    bragg_angle_corrected = diffraction_setup.angleBraggCorrected(energy)
    print("Bragg angle corrected for E=%f eV is %f deg" % (energy, numpy.degrees(bragg_angle_corrected)))

    darwin_half_width = diffraction_setup.darwinHalfwidth(energy)
    print("Darwin half width for E=%f eV is %f deg" % (energy, numpy.degrees(darwin_half_width[0])))


    deviations = numpy.linspace(angle_deviation_min, angle_deviation_max, angle_deviation_points)

    bunch_in = ComplexAmplitudePhotonBunch()
    for ia, deviation in enumerate(deviations):

        photon = ComplexAmplitudePhoton(energy_in_ev=energy,
                                        direction_vector=diffraction_setup.vectorIncomingPhotonDirection(
                                            energy,
                                            deviation * darwin_half_width[0],
                                            angle_center_flag=1,
                                        ),
                                        Esigma=1,
                                        Epi=1)

        bunch_in.addPhoton(photon)

    bunch_out = Diffraction.calculateDiffractedComplexAmplitudePhotonBunch(diffraction_setup,
                                                                           bunch_in,
                                                                           calculation_method=calculation_method,
                                                                           is_thick=is_thick,
                                                                           use_transfer_matrix=use_transfer_matrix,
                                                                           calculation_strategy_flag=calculation_strategy_flag,
                                                                           )

    bunch_out_dict = bunch_out.toDictionary()

    if do_plot:
        from srxraylib.plot.gol import plot
        plot(deviations,bunch_out_dict["intensity"],
             deviations,bunch_out_dict["intensityS"],
             deviations,bunch_out_dict["intensityP"],
             xtitle="y_zac",title="Reflectivity calculation using ComplexAmplitudePhoton calculation_method:%d"%calculation_method,
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
        is_thick=0,
        use_transfer_matrix=0,
        do_plot=0,
        calculation_strategy_flag=0,  # 0=mpmath 1=numpy 2=numpy-truncated
        dabax=None, # only used if material_constants_library_flag > 0
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
    print("is_thick                           = ", is_thick)
    print("use_transfer_matrix                = ", use_transfer_matrix)
    print("geometry_type_index                = ", geometry_type_index)
    if dabax is None:
        print("dabax                              = ", dabax)
    else:
        print("dabax                              = \n", dabax.info())

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
        dabax=dabax,
        )

    if theta is None:
        energy_mean = 0.5 * (energy_max + energy_min)
        theta = diffraction_setup.angleBragg(energy_mean)
        print("Using theta as the Bragg angle for mean E=%f eV, which is %f deg"%(energy_mean, numpy.degrees(theta)))

    energies = numpy.linspace(energy_min, energy_max, energy_points)

    bunch_in = ComplexAmplitudePhotonBunch()
    for ia, energy in enumerate(energies):
        bragg_angle = diffraction_setup.angleBragg(energy)
        deviation = theta - bragg_angle
        photon = ComplexAmplitudePhoton(energy_in_ev=energy,
                                        direction_vector=diffraction_setup.vectorIncomingPhotonDirection(
                                            energy,
                                            deviation,
                                            angle_center_flag=2,
                                        ),
                                        Esigma=1,
                                        Epi=1)

        bunch_in.addPhoton(photon)

    bunch_out = Diffraction.calculateDiffractedComplexAmplitudePhotonBunch(diffraction_setup,
                                                                           bunch_in,
                                                                           calculation_method=calculation_method,
                                                                           is_thick=is_thick,
                                                                           use_transfer_matrix=use_transfer_matrix,
                                                                           calculation_strategy_flag=calculation_strategy_flag)

    bunch_out_dict = bunch_out.toDictionary()

    if do_plot:
        from srxraylib.plot.gol import plot
        plot(energies,bunch_out_dict["intensity"],
             energies,bunch_out_dict["intensityS"],
             energies,bunch_out_dict["intensityP"],
             xtitle="photon energy [eV]",title="Reflectivity calculation using ComplexAmplitudePhoton calculation_method:%d"%calculation_method,
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
        is_thick=0,
        use_transfer_matrix=0,
        geometry_type_index=0,
        do_plot=0,
        calculation_strategy_flag=0,  # 0=mpmath 1=numpy 2=numpy-truncated
        dabax=None,  # only used if material_constants_library_flag > 0
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
    print("is_thick                           = ", is_thick)
    print("use_transfer_matrix                = ", use_transfer_matrix)
    print("geometry_type_index                = ", geometry_type_index)
    if dabax is None:
        print("dabax                              = ", dabax)
    else:
        print("dabax                              = \n", dabax.info())

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
        dabax=dabax,
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

    deviations = numpy.linspace(angle_deviation_min, angle_deviation_max, angle_deviation_points)
    energies   = numpy.linspace(energy_min, energy_max, energy_points)

    bunch_in = ComplexAmplitudePhotonBunch()
    print("Creating bunch...")
    t0 = time.time()
    for ie, energy in enumerate(energies):
        for ia, deviation in enumerate(deviations):
            photon = ComplexAmplitudePhoton(energy_in_ev=energy,
                                            direction_vector=diffraction_setup.vectorIncomingPhotonDirection(
                                                energy_mean,
                                                deviation,
                                                angle_center_flag=angle_center_flag,
                                            ),
                                            Esigma=1,
                                            Epi=1)

            bunch_in.addPhoton(photon)
    print("Done (%d ms)." % (1e3*(time.time() - t0)))

    print("Calculating diffraction...")
    t0 = time.time()
    bunch_out = Diffraction.calculateDiffractedComplexAmplitudePhotonBunch(diffraction_setup,
                                                                           bunch_in,
                                                                           calculation_method=calculation_method,
                                                                           is_thick=is_thick,
                                                                           use_transfer_matrix=use_transfer_matrix,
                                                                           calculation_strategy_flag=calculation_strategy_flag)

    print("Done (%d ms)." % (1e3*(time.time() - t0)))

    bunch_out_dict = bunch_out.toDictionary()

    if do_plot:
        from srxraylib.plot.gol import plot

        if angle_deviation_points == 1:
            plot(energies,bunch_out_dict["intensity"],
                 energies,bunch_out_dict["intensityS"],
                 energies,bunch_out_dict["intensityP"],
                 xtitle="photon energy [eV]",title="Reflectivity calculation using ComplexAmplitudePhoton calculation_method:%d"%calculation_method,
                 legend=["Total","Sigma","Pi"])

        if energy_points == 1:
            plot(1e6*deviations,bunch_out_dict["intensity"],
                 1e6*deviations,bunch_out_dict["intensityS"],
                 1e6*deviations,bunch_out_dict["intensityP"],
                 xtitle="theta - thetaB [urad]",title="Reflectivity calculation using ComplexAmplitudePhoton calculation_method:%d"%calculation_method,
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

    calculation_method = 1
    is_thick = 0
    use_transfer_matrix = 0
    calculation_strategy_flag = 0  # 0=mpmath 1=numpy 2=numpy-truncated

    if True:
        calc_xcrystal_angular_scan(material_constants_library_flag=0,
                                   do_plot=True,
                                   calculation_method=calculation_method,
                                   is_thick=is_thick,
                                   use_transfer_matrix=use_transfer_matrix,
                                   calculation_strategy_flag=calculation_strategy_flag, # 0=mpmath 1=numpy 2=numpy-truncated
                                   )

    if True:
        calc_xcrystal_angular_scan(material_constants_library_flag=0,
                                   geometry_type_index=1,
                                   thickness=10e-6,
                                   asymmetry_angle=numpy.radians(90),
                                   do_plot=True,
                                   calculation_method=calculation_method,
                                   is_thick=is_thick,
                                   use_transfer_matrix=use_transfer_matrix,
                                   calculation_strategy_flag=calculation_strategy_flag,
                                   )

        calc_xcrystal_energy_scan(material_constants_library_flag=0,
                                  do_plot=True,
                                  calculation_method=calculation_method,
                                  is_thick=is_thick,
                                  use_transfer_matrix=use_transfer_matrix,
                                  calculation_strategy_flag=calculation_strategy_flag,
                                  )

        calc_xcrystal_alphazachariasen_scan(do_plot=1,
                                            calculation_method=calculation_method,
                                            is_thick=is_thick,
                                            use_transfer_matrix=use_transfer_matrix,
                                            calculation_strategy_flag=calculation_strategy_flag,
                                            )

    if True:
        calc_xcrystal_double_scan(
            material_constants_library_flag=0,
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
            calculation_method=calculation_method,
            is_thick=is_thick,
            use_transfer_matrix=use_transfer_matrix,
            geometry_type_index=0,
            do_plot=1,
            calculation_strategy_flag=calculation_strategy_flag,
        )

        calc_xcrystal_double_scan(
            material_constants_library_flag=0,
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
            calculation_method=calculation_method,
            is_thick=is_thick,
            use_transfer_matrix=use_transfer_matrix,
            geometry_type_index=0,
            do_plot=1,
            calculation_strategy_flag=calculation_strategy_flag,
        )

    if True:
        calc_xcrystal_double_scan(
            material_constants_library_flag=0,
            crystal_name="Si",
            thickness=0.010,
            miller_h=1,
            miller_k=1,
            miller_l=1,
            asymmetry_angle=0.0,
            energy_min=7990,
            energy_max=8010,
            energy_points=150,
            angle_deviation_min=-100e-6,
            angle_deviation_max=100e-6,
            angle_deviation_points=150,
            angle_center_flag=2, # 0=Absolute angle, 1=Theta Bragg Corrected, 2=Theta Bragg
            calculation_method=calculation_method,
            is_thick=is_thick,
            use_transfer_matrix=use_transfer_matrix,
            geometry_type_index=0,
            do_plot=1,
            calculation_strategy_flag=calculation_strategy_flag,
            )