
#
# Study the validity of the Scherrer equation as compared with the dynamical diffraction
# Reproduces Fig. 5 in http://dx.doi.org/10.1107/S205327331600365X
#
import numpy


def calcFWHM(h, binSize):
    t = numpy.where(h >= max(h) * 0.5)
    return binSize * (t[0][-1] - t[0][0] + 1), t[0][-1], t[0][0]

def compute_crystalpy_profile(do_plot=0,
                              miller_h=1,
                              miller_k=1,
                              miller_l=1,
                              energy=8047.7969,
                              thickness=1e-08,
                              ):
    import numpy
    from crystalpy.util.calc_xcrystal import calc_xcrystal_angular_scan

    bunch_out_dict, diffraction_setup, deviations = calc_xcrystal_angular_scan(
        # material_constants_library_flag=self.material_constants_library_flag,
        crystal_name='LaB6_NIST',
        thickness=thickness,
        miller_h=miller_h,
        miller_k=miller_k,
        miller_l=miller_l,
        asymmetry_angle=0.0,
        energy=energy,
        angle_deviation_min=-0.001, # =-0.0523599,
        angle_deviation_max=0.001, # =0.0523599,
        angle_deviation_points=1000,
        angle_center_flag=2,
        calculation_method=0,
        is_thick=0,
        use_transfer_matrix=0,
        geometry_type_index=0,
        calculation_strategy_flag=0,
    )

    tmp = numpy.zeros((bunch_out_dict["energies"].size, 7))
    tmp[:, 0] = deviations / 0.017453292519943295
    tmp[:, 1] = energy
    tmp[:, 2] = bunch_out_dict["phaseP"]
    tmp[:, 3] = bunch_out_dict["phaseS"]
    # tmp[:, 4] = circular polarization
    tmp[:, 5] = bunch_out_dict["intensityP"]
    tmp[:, 6] = bunch_out_dict["intensityS"]

    if do_plot:
        from srxraylib.plot.gol import plot
        plot(tmp[:, 0], tmp[:, 6], tmp[:, 0], tmp[:, 5], xtitle="angle", legend=["S-pol", "P-pol"], title="%g [nm]"% (thickness*1e9))

    ba = numpy.degrees(diffraction_setup.angleBragg(energy))
    x = tmp[:, 0]
    y = tmp[:, 6]

    fwhm = calcFWHM(y,x[1]-x[0])
    fwhm = 2 * fwhm[0]

    return 2*ba, fwhm

if __name__ == "__main__":
    import  scipy.constants as codata

    # xtal_size = numpy.linspace(100, 10000, 10)
    xtal_size = numpy.array([200,400,500,600,700,800,900,1000,2000,5000,8000,10000])
    wavelength = 1.540598e-10
    photon_energy = codata.h * codata.c / codata.e / wavelength
    print("photon_energy is: ", photon_energy)

    RA110 = numpy.zeros_like(xtal_size, dtype=float)
    RA333 = numpy.zeros_like(xtal_size, dtype=float)
    for i in range(xtal_size.size):
        try:
            ba, fwhm = compute_crystalpy_profile(do_plot=0,miller_h=1, miller_k=1, miller_l=0,
                                                 thickness=xtal_size[i]*1e-9)
            theta_bragg = numpy.radians(ba / 2)
            D = xtal_size[i]*1e-9
            scherrer_fwhm = 4 * numpy.sqrt(numpy.log(2) / numpy.pi) * wavelength * numpy.abs(
                numpy.sin(theta_bragg)) / D / numpy.sin(2 * theta_bragg)
            print('scherrer, fwhm, ratio: ', numpy.degrees(scherrer_fwhm), fwhm, numpy.degrees(scherrer_fwhm) / fwhm)
            RA110[i] = numpy.degrees(scherrer_fwhm) / fwhm
        except:
            pass

    for i in range(xtal_size.size):
        try:
            ba, fwhm = compute_crystalpy_profile(do_plot=0,miller_h=3, miller_k=3, miller_l=3,
                                                 thickness=xtal_size[i]*1e-9)
            theta_bragg = numpy.radians(ba / 2)
            D = xtal_size[i]*1e-9
            scherrer_fwhm = 4 * numpy.sqrt(numpy.log(2) / numpy.pi) * wavelength * numpy.abs(
                numpy.sin(theta_bragg)) / D / numpy.sin(2 * theta_bragg)
            print('scherrer, fwhm, ratio: ', numpy.degrees(scherrer_fwhm), fwhm, numpy.degrees(scherrer_fwhm) / fwhm)
            RA333[i] = numpy.degrees(scherrer_fwhm) / fwhm
        except:
            pass


    from srxraylib.plot.gol import plot
    plot(xtal_size, RA110,
         xtal_size, RA333,
         xlog=1, marker=['+','x'], legend=['LaB6 110', 'LaB6 333'], yrange=[0,1.2], xrange=[10,10000],
         xtitle="crystal size [nm]", ytitle="FWHM_Scherrer/FWHM_DynDiff",
         # BA, FWHM, marker=[None,'+'], linestyle=[None,''], xtitle="2 theta [deg]", ytitle="FWHM [deg]",
         # legend=['Scherrer eq. (3) in Leitao Muniz et al.', 'dynamical diffraction with crystalpy']
         )

