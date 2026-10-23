"""
Unittest for the Stokes parameters calculated by calc_xcrystal_angular_scan (flag_calculate_stokes > 0), and their
convolution with a Gaussian (flag_calculate_stokes_convolution=1).

Conventions of "Producing and Measuring Circularly Polarized X-rays: A Jones and Mueller-Stokes Treatment of Crystal
Phase Plates": time factor exp(-i omega t), J(chi) = R(-chi) diag(rs, rp) R(chi), S3 = 2 Im(Ex* Ey) > 0 for LCP.
"""
import unittest
import contextlib, io

import numpy
from numpy.testing import assert_allclose

from crystalpy.util.calc_xcrystal import calc_xcrystal_angular_scan, apply_convolution_with_gaussian

# Si 111 phase plate (transmitted beam in Bragg geometry) at 2838 eV, 10 um, rotated by chi=-45 deg
FWHM = 20e-6 # rad


def scan(calculation_method=1, flag_calculate_stokes=1, convolution=0, fwhm=0.0, jones_in=(1, 0), chi_deg=-45.0,
         npoints=601):
    with contextlib.redirect_stdout(io.StringIO()):
        bunch_out_dict, diffraction_setup, deviations = calc_xcrystal_angular_scan(
            crystal_name="Si", thickness=10e-6, miller_h=1, miller_k=1, miller_l=1, asymmetry_angle=0.0,
            energy=2838.0, angle_deviation_min=-600e-6, angle_deviation_max=600e-6,
            angle_deviation_points=npoints, angle_center_flag=1,
            flag_calculate_stokes=flag_calculate_stokes, jones_in=list(jones_in), chi_deg=chi_deg,
            flag_calculate_stokes_convolution=convolution, fwhm_calculate_stokes_convolution=fwhm,
            calculation_method=calculation_method, is_thick=0, use_transfer_matrix=0,
            geometry_type_index=2, calculation_strategy_flag=0, material_constants_library_flag=0)
    return bunch_out_dict, deviations


def stokes(d):
    return numpy.array([numpy.asarray(d[k]) for k in ("P0", "P1", "P2", "P3")])


class CalcXcrystalStokesTest(unittest.TestCase):

    def test_jones_and_stokes(self):
        # Stokes parameters from the Jones vector J(chi) E_in, with J(chi) = R(-chi) diag(rs, rp) R(chi)
        E_in = numpy.array([0.6, 0.8 * numpy.exp(0.7j)])
        for chi_deg in [-45.0, 30.0]:
            d, _ = scan(jones_in=E_in, chi_deg=chi_deg, npoints=41)
            rs, rp = d["complexAmplitudeS"], d["complexAmplitudeP"]
            c, s = numpy.cos(numpy.radians(chi_deg)), numpy.sin(numpy.radians(chi_deg))
            R = lambda c, s: numpy.array([[c, s], [-s, c]])
            Ex = numpy.zeros(rs.size, dtype=complex)
            Ey = numpy.zeros(rs.size, dtype=complex)
            for i in range(rs.size):
                E = R(c, -s) @ numpy.diag([rs[i], rp[i]]) @ R(c, s) @ E_in
                Ex[i], Ey[i] = E
            expected = numpy.array([abs(Ex) ** 2 + abs(Ey) ** 2, abs(Ex) ** 2 - abs(Ey) ** 2,
                                    2 * numpy.real(Ex * numpy.conj(Ey)), 2 * numpy.imag(numpy.conj(Ex) * Ey)])
            assert_allclose(stokes(d), expected, atol=1e-14)

    def test_zachariasen_conjugated(self):
        # the Zachariasen amplitudes (exp(+i omega t)) are conjugated: same Stokes parameters as Guigay
        g, _ = scan(calculation_method=1, npoints=81)
        z, _ = scan(calculation_method=0, npoints=81)
        assert_allclose(stokes(z), stokes(g), atol=2e-3)
        self.assertGreater(numpy.abs(stokes(g)[3]).max(), 0.01) # S3 is not negligible (the test is meaningful)

    def test_convolution(self):
        d0, deviations = scan()
        S = stokes(d0)
        sigma = FWHM / (2 * numpy.sqrt(2 * numpy.log(2)))
        S_conv = numpy.array([apply_convolution_with_gaussian(deviations, Sk, sigma=sigma) for Sk in S])

        # non normalized Stokes parameters
        d1, _ = scan(convolution=1, fwhm=FWHM)
        assert_allclose(stokes(d1), S_conv, atol=1e-15)
        self.assertEqual(d1["stokes_convolution_fwhm"], FWHM)

        # normalized: the non normalized ones are convolved, and then normalized
        d2, _ = scan(flag_calculate_stokes=2, convolution=1, fwhm=FWHM)
        assert_allclose(stokes(d2)[0], S_conv[0], atol=1e-15)
        assert_allclose(stokes(d2)[1:], S_conv[1:] / S_conv[0], atol=1e-12)

        # the beam becomes partially polarized close to the Bragg angle
        P = stokes(d2)[1:]
        DOP = numpy.sqrt(numpy.sum(P ** 2, axis=0))
        self.assertLess(DOP.min(), 0.9)
        assert_allclose(DOP[[0, -1]], 1.0, atol=1e-3) # far from the Bragg angle

        # no artifacts at the edges (the curves are extended with their end values)
        assert_allclose(stokes(d1)[0][[0, -1]], S[0][[0, -1]], rtol=1e-2)

        # no convolution for fwhm=0
        d3, _ = scan(convolution=1, fwhm=0.0)
        assert_allclose(stokes(d3), S, atol=0)

    def test_gaussian_normalization(self):
        x = numpy.linspace(-1, 1, 2001)
        y = numpy.ones_like(x)
        assert_allclose(apply_convolution_with_gaussian(x, y, sigma=0.05), y, atol=1e-12)
        y = numpy.where(numpy.abs(x) < 0.001, 1.0 / (x[1] - x[0]), 0.0) # discrete delta (area 1)
        yc = apply_convolution_with_gaussian(x, y, sigma=0.05)
        assert_allclose(numpy.sum(yc) * (x[1] - x[0]), 1.0, atol=1e-6)
        assert_allclose(numpy.sqrt(numpy.sum(yc * x ** 2) / numpy.sum(yc)), 0.05, rtol=1e-3)


if __name__ == "__main__":
    unittest.main()
