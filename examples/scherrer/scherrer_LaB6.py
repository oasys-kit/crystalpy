#
# Study the validity of the Scherrer equation as compared with the dynamical diffraction
# Reproduces Fig. 1 in http://dx.doi.org/10.1107/S205327331600365X
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
                              crystal_name='LaB6_NIST',
                              ):
    import numpy
    from crystalpy.util.calc_xcrystal import calc_xcrystal_angular_scan, calc_xcrystal_energy_scan, \
        calc_xcrystal_alphazachariasen_scan

    bunch_out_dict, diffraction_setup, deviations = calc_xcrystal_angular_scan(
        # material_constants_library_flag=self.material_constants_library_flag,
        crystal_name=crystal_name,
        thickness=1e-08,
        miller_h=miller_h,
        miller_k=miller_k,
        miller_l=miller_l,
        asymmetry_angle=0.0,
        energy=energy,
        angle_deviation_min=-0.11, # =-0.0523599,
        angle_deviation_max=0.11, # =0.0523599,
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
        plot(tmp[:, 0], tmp[:, 6], tmp[:, 0], tmp[:, 5], xtitle="angle", legend=["S-pol", "P-pol"], title="%d%d%d"% (miller_h,miller_k,miller_l))

    ba = numpy.degrees(diffraction_setup.angleBragg(energy))
    x = tmp[:, 0]
    y = tmp[:, 6]

    fwhm = calcFWHM(y,x[1]-x[0])
    fwhm = 2 * fwhm[0]

    return 2*ba, fwhm

if __name__ == "__main__":
    import  scipy.constants as codata

    theta_bragg_deg = numpy.linspace(10, 70, 100)
    theta_bragg = numpy.radians( theta_bragg_deg )
    D = 10e-9 # 10 nm
    # D = 1e-6 # 1 um
    wavelength = 1.540598e-10
    photon_energy = codata.h * codata.c / codata.e / wavelength
    print("photon_energy is: ", photon_energy)

    scherrer_fwhm = 4 * numpy.sqrt(numpy.log(2) / numpy.pi) * wavelength * numpy.abs(numpy.sin(theta_bragg)) / D / numpy.sin(2 * theta_bragg)

    # print("D: %g, wavelength:%g, FWHM:%g" % (D, wavelength, numpy.degrees(scherrer_fwhm)))
    #

    H = [1,1,1,2,2,2,2,2,2,3,3,2,3,3,4,4,4,3,4,4,3,4,5,5,5]
    K = [0,1,1,0,0,1,1,2,2,1,1,2,2,2,0,1,1,3,2,2,3,2,0,1,1]
    L = [0,0,1,0,0,0,1,0,1,0,1,2,0,1,0,0,1,1,0,1,2,2,0,0,1]

    BA  = numpy.zeros(len(H))
    FWHM = numpy.zeros(len(H))
    FWHM_S = numpy.zeros(len(H))

    crystal_name = 'LaB6_NIST'


    for i in range(len(H)):
        try:
            ba, fwhm = compute_crystalpy_profile(do_plot=0,miller_h=H[i], miller_k=K[i], miller_l=L[i],
                                                 crystal_name=crystal_name)
            BA[i] = ba
            FWHM[i] = fwhm
            theta_bragg1 = numpy.radians(ba / 2)
            print(">>>>>2theta: ", ba, 'FWHM: ', fwhm, 'thetaB: ', numpy.degrees(theta_bragg1))
            FWHM_S[i] = 4 * numpy.sqrt(numpy.log(2) / numpy.pi) * wavelength * numpy.abs(numpy.sin(theta_bragg1)) / D / numpy.sin(2 * theta_bragg1)
        except:
            print(">>>> ERROR FOUND! h=%d,k=%d,l=%d" % (H[i], K[i], L[i]))
            pass

    # print(ba,fwhm)

    from srxraylib.plot.gol import plot
    plot(2 * theta_bragg_deg, numpy.degrees(scherrer_fwhm),
         BA, FWHM, marker=[None,'+'], linestyle=[None,''], xtitle="2 theta [deg]", ytitle="FWHM [deg]",
         legend=['Scherrer eq. (3) in Leitao Muniz et al.', 'dynamical diffraction with crystalpy'],
         title=crystal_name,
         )

    print(BA, FWHM, FWHM_S)
