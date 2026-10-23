"""
Utility functions for crystal diffraction calculations (xcrystal tools).
"""
import numpy
import numpy as np
import time
import copy

from scipy.signal import fftconvolve

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
        angle_center_flag=2,     # 0=Absolute angle, 1=Theta Bragg Corrected, 2=Theta Bragg
        flag_calculate_stokes=0, # **new** 0=No Stokes calculation, 1=Yes (S0, S1, S2, S3), 2=Yes, normalized (P1=S1/S0...)
        chi_deg=0.0,     # for flag_calculate_stokes>0: rotation of the crystal around the beam (from x toward y)
        jones_in=[1,1], # for flag_calculate_stokes>0: Jones vector (Ex, Ey) of the incident beam
        flag_calculate_stokes_convolution=0, # for flag_calculate_stokes>0: 1=convolve the Stokes parameters with a Gaussian
        fwhm_calculate_stokes_convolution=0.0, # for flag_calculate_stokes_convolution=1: FWHM of the Gaussian (rad)
        calculation_method=0,
        is_thick=0,
        use_transfer_matrix=0,
        do_plot=0,
        calculation_strategy_flag=0,  # 0=mpmath 1=numpy 2=numpy-truncated
        dabax=None, # only used if material_constants_library_flag > 0
        ):
    """
    Angular scan of the diffraction (or transmission) of a crystal for a monochromatic plane wave, and (optionally)
    the Stokes parameters of the diffracted (or transmitted) beam.

    Parameters
    ----------
    crystal_name : str, optional
        The crystal name (in the xraylib or DABAX crystal list).
    miller_h, miller_k, miller_l : int, optional
        The Miller indices.
    thickness : float, optional
        The crystal thickness in m.
    asymmetry_angle : float, optional
        The asymmetry angle in rad.
    geometry_type_index : int, optional
        0=Bragg diffraction, 1=Laue diffraction, 2=Bragg transmission, 3=Laue transmission.
    material_constants_library_flag : int, optional
        0=xraylib, 1=DABAX (with the dabax instance).
    energy : float, optional
        The photon energy in eV.
    angle_deviation_min, angle_deviation_max : float, optional
        The limits of the angular scan in rad (deviations from the angle defined by angle_center_flag).
    angle_deviation_points : int, optional
        The number of points of the scan.
    angle_center_flag : int, optional
        0=absolute angle, 1=deviation from the Bragg angle corrected for refraction, 2=deviation from the Bragg angle.
    flag_calculate_stokes : int, optional
        0=No; 1=Stokes parameters S0, S1, S2, S3 of the output beam (stored in the keys "P0"..."P3");
        2=S0 and the normalized P1=S1/S0, P2=S2/S0, P3=S3/S0. The conventions are those of "Producing and Measuring
        Circularly Polarized X-rays: A Jones and Mueller-Stokes Treatment of Crystal Phase Plates": time factor
        exp(-i omega t), S3 = 2 Im(Ex* Ey) > 0 for left-handed polarization, J(chi) = R(-chi) diag(rs, rp) R(chi)
        with R(chi) = [[cos, sin], [-sin, cos]].
    chi_deg : float, optional
        For flag_calculate_stokes > 0, the rotation of the crystal around the beam in deg (positive from x toward y).
        Note: the orientation angle alpha of a SHADOW4 element is alpha = -chi.
    jones_in : list, optional
        For flag_calculate_stokes > 0, the Jones vector (Ex, Ey) of the incident beam, e.g., [1, 0] for horizontal
        linear polarization.
    flag_calculate_stokes_convolution : int, optional
        For flag_calculate_stokes > 0: 1=convolve the Stokes parameters versus the angular deviation with a
        normalized Gaussian, to approximate an incident beam made of incoherent plane waves with a Gaussian
        distribution of deviations from the Bragg angle (e.g., from its divergence and, through
        (Delta E / E) tan(theta_B), its energy spread). The non-normalized S0...S3 are convolved (they are additive
        for incoherent waves), and then normalized if flag_calculate_stokes=2. The curves are extended with their end
        values at the edges of the scan. The scan must be equally spaced and wider than the Gaussian.
    fwhm_calculate_stokes_convolution : float, optional
        For flag_calculate_stokes_convolution=1, the FWHM of the Gaussian in rad (sigma = FWHM / (2 sqrt(2 ln 2))).
    calculation_method : int, optional
        0=Zachariasen, 1=Guigay. The Zachariasen amplitudes (time factor exp(+i omega t)) are conjugated for the
        calculation of the Stokes parameters; the amplitudes and phases stored in the output dictionary are not.
    is_thick : int, optional
        1=use the thick crystal approximation.
    use_transfer_matrix : int, optional
        1=use the transfer matrix (Guigay method).
    do_plot : int, optional
        1=plot the result.
    calculation_strategy_flag : int, optional
        0=mpmath, 1=numpy, 2=numpy-truncated.
    dabax : instance of DabaxXraylib, optional
        The DABAX access (for material_constants_library_flag=1).

    Returns
    -------
    tuple
        (bunch_out_dict, diffraction_setup, deviations): the dictionary with the results (keys "intensityS",
        "intensityP", "phaseS", "phaseP", "complexAmplitudeS", "complexAmplitudeP", "deviations", and, if
        flag_calculate_stokes > 0, "P0", "P1", "P2", "P3", "stokes_normalized", "stokes_convolution_fwhm"), the
        diffraction setup, and the
        array of angular deviations (rad).
    """

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
    print("chi_deg                            = ", chi_deg)
    print("flag_calculate_stokes              = ", flag_calculate_stokes)
    print("chi_deg                            = ", chi_deg)
    print("jones_in                           = ", jones_in)
    print("flag_calculate_stokes_convolution  = ", flag_calculate_stokes_convolution)
    print("fwhm_calculate_stokes_convolution  = ", fwhm_calculate_stokes_convolution)
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

    if flag_calculate_stokes > 0:
        #
        # Stokes parameters of the diffracted (or transmitted) beam for an incident beam with the Jones vector
        # jones_in = (Ex, Ey), and the crystal rotated around the beam by chi_deg.
        #
        # Conventions of "Producing and Measuring Circularly Polarized X-rays: A Jones and Mueller-Stokes Treatment
        # of Crystal Phase Plates" (paper-polarisation, main11):
        #   * time factor exp(-i omega t) [physical field Re((Ex x + Ey y) exp(-i(omega t - k z)))], right-handed
        #     basis (x, y, z=k). The amplitudes of the Guigay method (calculation_method=1) are in this convention;
        #     those of the Zachariasen method (calculation_method=0) use exp(+i omega t), and they are complex
        #     conjugated here (conjugation changes the sign of S3 and of the retardation phi_p - phi_s).
        #   * crystal Jones matrix J0 = diag(rs, rp), with sigma along x and pi along y for chi=0;
        #     retardation Delta phi = arg(rp) - arg(rs).
        #   * rotation of the crystal by chi around the beam: J(chi) = R(-chi) J0 R(chi), with the passive
        #     rotation R(chi) = [[cos, sin], [-sin, cos]]; positive chi rotates the crystal sigma axis from x toward y.
        #     The output Jones vector is in the same (laboratory) basis as the input.
        #   * Stokes parameters: S0 = |Ex|^2 + |Ey|^2, S1 = |Ex|^2 - |Ey|^2, S2 = 2 Re(Ex Ey*), S3 = 2 Im(Ex* Ey)
        #     (i.e., Sk = Tr(sigma_k rho)); S3 > 0 is left-handed (LCP), e.g., (1, +i)/sqrt(2) gives S3 = +1 (LCP)
        #     and (1, -i)/sqrt(2) gives S3 = -1 (RCP).
        # The results are stored in the keys "P0", "P1", "P2", "P3": the Stokes parameters S0...S3 for
        # flag_calculate_stokes=1, or S0 and the normalized P1=S1/S0, P2=S2/S0, P3=S3/S0 for flag_calculate_stokes=2.
        #
        # Note (SHADOW4): the same conventions are used by SHADOW4 (time factor, standard basis (x, y, k) =
        # (u_sigma, -u_pi, v) with u_pi = u_sigma x v, Stokes parameters), but the orientation angle alpha of a SHADOW4
        # element is measured from u_sigma toward u_pi, i.e., in the opposite sense: alpha = -chi.
        #
        chi = numpy.radians(chi_deg)

        rs = bunch_out_dict["complexAmplitudeS"]
        rp = bunch_out_dict["complexAmplitudeP"]
        if calculation_method == 0: # Zachariasen: exp(+i omega t) -> exp(-i omega t)
            rs = numpy.conjugate(rs)
            rp = numpy.conjugate(rp)

        sa = numpy.sin(chi)
        ca = numpy.cos(chi)
        cs = sa * ca
        c2 = ca ** 2
        s2 = sa ** 2

        # [[r00, r01], [r10, r11]] = R(-chi) [[rs,0],[0,rp]] R(chi)
        r00 = rs * c2 + rp * s2
        r01 = cs * (rs - rp)
        r10 = cs * (rs - rp)
        r11 = rs * s2 + rp * c2
        jones_out = (r00 * jones_in[0] + r01 * jones_in[1], r10 * jones_in[0] + r11 * jones_in[1])

        P0_out = numpy.abs(jones_out[0]) ** 2 + numpy.abs(jones_out[1]) ** 2
        P1_out = numpy.abs(jones_out[0]) ** 2 - numpy.abs(jones_out[1]) ** 2
        P2_out = 2 * numpy.real( jones_out[0] * numpy.conjugate(jones_out[1]) )
        P3_out = 2 * numpy.imag( numpy.conjugate(jones_out[0]) * jones_out[1] )

        # convolution of the (non normalized) Stokes parameters with a Gaussian of FWHM fwhm_calculate_stokes_convolution
        if flag_calculate_stokes_convolution and fwhm_calculate_stokes_convolution > 0:
            if angle_deviation_points < 2:
                raise ValueError("The convolution of the Stokes parameters needs an angular scan (angle_deviation_points > 1)")
            sigma = fwhm_calculate_stokes_convolution / (2 * numpy.sqrt(2 * numpy.log(2)))
            P0_out = apply_convolution_with_gaussian(deviations, P0_out, sigma=sigma)
            P1_out = apply_convolution_with_gaussian(deviations, P1_out, sigma=sigma)
            P2_out = apply_convolution_with_gaussian(deviations, P2_out, sigma=sigma)
            P3_out = apply_convolution_with_gaussian(deviations, P3_out, sigma=sigma)
            bunch_out_dict["stokes_convolution_fwhm"] = fwhm_calculate_stokes_convolution
        else:
            bunch_out_dict["stokes_convolution_fwhm"] = 0.0

        if flag_calculate_stokes == 2:
            P1_out = P1_out / P0_out
            P2_out = P2_out / P0_out
            P3_out = P3_out / P0_out

        bunch_out_dict["P0"] = P0_out
        bunch_out_dict["P1"] = P1_out
        bunch_out_dict["P2"] = P2_out
        bunch_out_dict["P3"] = P3_out
        bunch_out_dict["stokes_normalized"] = (flag_calculate_stokes == 2)

    bunch_out_dict["deviations"] = deviations

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

def apply_convolution_with_gaussian(x, y, sigma=7.0):
    """
    Convolution of a curve y(x) (x equally spaced) with a normalized Gaussian of standard deviation sigma (same
    units as x). The curve is extended with its end values to avoid the artifacts of the zero padding at the edges.

    Parameters
    ----------
    x : numpy array
        The abscissas (equally spaced).
    y : numpy array
        The ordinates.
    sigma : float
        The rms width of the Gaussian, in the units of x. If zero, y is returned.

    Returns
    -------
    numpy array
        The convolved curve, at the same abscissas.
    """
    if sigma <= 0: return numpy.array(y, copy=True)
    dx = x[1] - x[0]
    npad = int(numpy.ceil(5 * sigma / abs(dx)))
    xk = numpy.arange(-npad, npad + 1) * dx
    gaussian = numpy.exp(-xk ** 2 / (2 * sigma ** 2))
    gaussian /= gaussian.sum()
    y_padded = numpy.concatenate((numpy.full(npad, y[0]), y, numpy.full(npad, y[-1])))
    return fftconvolve(y_padded, gaussian, mode="same")[npad:npad + len(y)]


def apply_convolution_on_bunch_out_dict(bunch_out_dict, sigma=0.1, deviations=None):
    """
    Convolution of the Stokes parameters of an angular scan (from calc_xcrystal_angular_scan with
    flag_calculate_stokes > 0) with a Gaussian of standard deviation sigma (in the units of the deviations, rad).

    The (non normalized) Stokes parameters S0...S3 are convolved, as they are additive for incoherent
    superposition; the normalized ones are then P1=S1/S0... (if the input was normalized).

    Parameters
    ----------
    bunch_out_dict : dict
        The dictionary returned by calc_xcrystal_angular_scan with flag_calculate_stokes > 0.
    sigma : float
        The rms width of the Gaussian (rad).
    deviations : numpy array, optional
        The angular deviations (rad); if None, bunch_out_dict["deviations"] is used.

    Returns
    -------
    dict
        A copy of bunch_out_dict with the convolved "P0", "P1", "P2", "P3".
    """
    bunch_out_dict_conv = copy.deepcopy(bunch_out_dict)
    if deviations is None: deviations = bunch_out_dict["deviations"]
    normalized = bunch_out_dict.get("stokes_normalized", False)
    S0 = numpy.asarray(bunch_out_dict["P0"])
    S = [numpy.asarray(bunch_out_dict[k]) for k in ("P1", "P2", "P3")]
    if normalized: S = [Sk * S0 for Sk in S]

    S0_conv = apply_convolution_with_gaussian(deviations, S0, sigma=sigma)
    S_conv = [apply_convolution_with_gaussian(deviations, Sk, sigma=sigma) for Sk in S]
    if normalized: S_conv = [Sk / S0_conv for Sk in S_conv]

    bunch_out_dict_conv["P0"] = S0_conv
    bunch_out_dict_conv["P1"], bunch_out_dict_conv["P2"], bunch_out_dict_conv["P3"] = S_conv
    return bunch_out_dict_conv

if __name__ == "__main__":
    from srxraylib.plot.gol import set_qt
    try:
        import xraylib
        material_constants_library_flag = 0
        dabax=None
    except:
        material_constants_library_flag = 1
        from dabax.dabax_xraylib import DabaxXraylib
        dabax = DabaxXraylib()

    calculation_method = 1
    is_thick = 0
    use_transfer_matrix = 0
    calculation_strategy_flag = 0  # 0=mpmath 1=numpy 2=numpy-truncated

    if True:
        calc_xcrystal_angular_scan(do_plot=True,
                                   calculation_method=calculation_method,
                                   is_thick=is_thick,
                                   use_transfer_matrix=use_transfer_matrix,
                                   calculation_strategy_flag=calculation_strategy_flag, # 0=mpmath 1=numpy 2=numpy-truncated
                                   material_constants_library_flag=material_constants_library_flag,
                                   dabax=dabax,
                                   )

    if True:
        calc_xcrystal_angular_scan(geometry_type_index=1,
                                   thickness=10e-6,
                                   asymmetry_angle=numpy.radians(90),
                                   do_plot=True,
                                   calculation_method=calculation_method,
                                   is_thick=is_thick,
                                   use_transfer_matrix=use_transfer_matrix,
                                   calculation_strategy_flag=calculation_strategy_flag,
                                   material_constants_library_flag=material_constants_library_flag,
                                   dabax=dabax,
                                   )

        calc_xcrystal_energy_scan(do_plot=True,
                                  calculation_method=calculation_method,
                                  is_thick=is_thick,
                                  use_transfer_matrix=use_transfer_matrix,
                                  calculation_strategy_flag=calculation_strategy_flag,
                                  material_constants_library_flag=material_constants_library_flag,
                                  dabax=dabax,
                                  )

        calc_xcrystal_alphazachariasen_scan(do_plot=1,
                                            calculation_method=calculation_method,
                                            is_thick=is_thick,
                                            use_transfer_matrix=use_transfer_matrix,
                                            calculation_strategy_flag=calculation_strategy_flag,
                                            material_constants_library_flag=material_constants_library_flag,
                                            dabax=dabax,
                                            )

    if True:
        calc_xcrystal_double_scan(
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
            material_constants_library_flag=material_constants_library_flag,
            dabax=dabax,
        )

        calc_xcrystal_double_scan(
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
            material_constants_library_flag=material_constants_library_flag,
            dabax=dabax,
        )

    if True:
        calc_xcrystal_double_scan(
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
            material_constants_library_flag=material_constants_library_flag,
            dabax=dabax,
            )