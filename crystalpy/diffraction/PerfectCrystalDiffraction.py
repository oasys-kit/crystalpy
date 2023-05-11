"""
Calculates crystal diffraction according to Zachariasen's representation of the dynamic theory of crystal diffraction
for perfect crystals.
Except for energy all units are in SI. Energy is in eV.
"""

import numpy
from crystalpy.util.Photon import Photon
from crystalpy.util.PhotonBunch import PhotonBunch
from crystalpy.diffraction.GeometryType import BraggDiffraction, LaueDiffraction, BraggTransmission, LaueTransmission

# Use mpmath if possible. Otherwise use native cmath.
try:
    # raise ImportError
    import mpmath
    use_mpmath = True
    mpmath_sin   =  numpy.vectorize(mpmath.sin)
    mpmath_cos   =  numpy.vectorize(mpmath.cos)
    mpmath_exp   =  numpy.vectorize(mpmath.exp)

    # mpmath_sin   =  (mpmath.sin)
    # mpmath_cos   =  (mpmath.cos)
    # mpmath_exp   =  (mpmath.exp)

except ImportError:
    # import cmath
    import numpy as cmath
    use_mpmath = False
    print("mpmath module for arbitrary-precision floating-point arithmetic could not be found!\n "
          "Using numpy instead. This could lead to overflow errors for Zachariasen calculations.\n")

if use_mpmath:
    print("Using mpmath.")
else:
    print("Using numpy instead of mpmath (not available).")


class CalculationStrategy(object):
    """
    Abstract strategy for calculation. Can be plain python or arbitrary precision like mpmath.
    """
    def createVariable(self, initial_value):
        """
        Factory method for calculation variable.
        :param initial_value: Initial value of the variable.
        :return: Calculation variable.
        """
        raise Exception("Must override this method.")

    def exponentiate(self, power):
        """
        Exponentiates to the power.
        :param power: The power to raise to.
        :return: Exponential.
        """
        raise Exception("Must override this method.")

    def sin(self, power):
        """
        Sin to the power.
        :param power: The power to raise to.
        :return: Sin.
        """
        raise Exception("Must override this method.")

    def cos(self, power):
        """
        Cos to the power.
        :param power: The power to raise to.
        :return: Cos.
        """
        raise Exception("Must override this method.")


    def toComplex(self, variable):
        """
        Converts calculation variable to native python complex.
        :param variable: Calculation variable to convert.
        :return: Native python complex variable.
        """
        raise Exception("Must override this method.")


class CalculationStrategyMPMath(CalculationStrategy):
    """
    Use mpmath for calculation.
    """
    def __init__(self):
        """
        Constructor.
        """
        # Use 32 digits in mpmath calculations.
        mpmath.mp.dps = 32

    def createVariable(self, initial_value):
        """
        Factory method for calculation variable.
        :param initial_value: Initial value of the variable.
        :return: mpmath variable.
        """

        if initial_value.size == 1:
            mpc = mpmath.mpc(complex(initial_value.real) + 1j * complex(initial_value.imag))
        else:
            mpc = mpmath.mpc(complex(1) + 1j * complex(0)) * initial_value

        return mpc

    def exponentiate(self, power):
        """
        Exponentiates to the power.
        :param power: The power to raise to.
        :return: Exponential.
        """
        # return mpmath.exp(power)
        if False: # isinstance(power, numpy.ndarray):
            power_out = numpy.zeros_like(power, dtype=mpmath.mpc)
            for i in range(power.size):
                power_out[i] = mpmath_exp(power[i])
            return power_out
        else:
            return mpmath_exp(power)

    def sin(self, power):
        """
        Sin to the power.
        :param power: The power to raise to.
        :return: Sin.
        """
        # return mpmath.sin(power)
        if False: # isinstance(power, numpy.ndarray):
            power_out = numpy.zeros_like(power, dtype=mpmath.mpc)
            for i in range(power.size):
                power_out[i] = mpmath_sin(power[i])
            return power_out
        else:
            return mpmath_sin(power)

    def cos(self, power):
        """
        Cos to the power.
        :param power: The power to raise to.
        :return: Cos.
        """
        # return mpmath.cos(power)
        if False: # isinstance(power, numpy.ndarray):
            power_out = numpy.zeros_like(power, dtype=mpmath.mpc)
            for i in range(power.size):
                power_out[i] = mpmath_cos(power[i])
            return power_out
        else:
            return mpmath_cos(power)

    def toComplex(self, variable):
        """
        Converts calculation variable to native python complex.
        :param variable: Calculation variable to convert.
        :return: Native python complex variable.
        """
        # return complex(variable)

        return numpy.array(variable, dtype=complex)


class CalculationStrategyMath(CalculationStrategy):
    """
    Use plain python for calculation.
    """
    def createVariable(self, initial_value):
        """
        Factory method for calculation variable.
        :param initial_value: Initial value of the variable.
        :return: mpmath variable.
        """
        return initial_value + 0j # complex(initial_value)

    def exponentiate(self, power):
        """
        Exponentiates to the power.
        :param power: The power to raise to.
        :return: Exponential.
        """
        try:
            ans =  cmath.exp(power)
        except:
            ans = float("Inf")
        return ans

    def sin(self, power):
        """
        Sin to the power.
        :param power: The power to raise to.
        :return: Sin.
        """
        return cmath.sin(power)

    def cos(self, power):
        """
        Cos to the power.
        :param power: The power to raise to.
        :return: Cos.
        """
        return cmath.cos(power)


    def toComplex(self, variable):
        """
        Converts calculation variable to native python complex.
        :param variable: Calculation variable to convert.
        :return: Native python complex variable.
        """
        return complex(variable)


class PerfectCrystalDiffraction(object):
    isDebug = False

    def __init__(self, geometry_type, bragg_normal, surface_normal, bragg_angle, psi_0, psi_H, psi_H_bar, thickness, d_spacing):
        """
        Constructor.
        :param geometry_type: The diffraction geometry, i.e. BraggDiffraction, LaueTransmission,...
        :param bragg_normal: Normal on the reflection planes.
        :param surface_normal:Norm on crystal surface pointing outward.
        :param bragg_angle: Bragg angle.
        :param psi_0: Psi0 as defined in Zachariasen [3-95].
        :param psi_H: PsiH as defined in Zachariasen [3-95].
        :param psi_H_bar: PsiHBar as defined in Zachariasen [3-95].
        :param thickness: Thickness of the crystal.
        :param d_spacing: Spacing of parallel planes.
        """
        self._geometryType = geometry_type
        self._bragg_normal = bragg_normal
        self._surface_normal = surface_normal
        self._bragg_angle = bragg_angle
        self._psi_0 = psi_0
        self._psi_H = psi_H
        self._psi_H_bar = psi_H_bar
        self._thickness = thickness
        self._d_spacing = d_spacing

        global use_mpmath
        if use_mpmath:
            self._calculation_strategy = CalculationStrategyMPMath()
        else:
            self._calculation_strategy = CalculationStrategyMath()

    def braggNormal(self):
        """
        Returns the Bragg normal, i.e. normal on the reflection planes.
        :return: Bragg normal.
        """
        return self._bragg_normal

    def surface_normal(self):
        """
        Returns the surface normal that points outwards the crystal.
        :return: Surface normal.
        """
        return self._surface_normal

    def braggAngle(self):
        """
        Returns the Bragg angle.
        :return: The Bragg angle.
        """
        return self._bragg_angle

    def Psi0(self):
        """
        Returns Psi0 as defined in Zachariasen [3-95].
        :return: Psi0.
        """
        return self._psi_0

    def PsiH(self):
        """
        Returns Psi0 as defined in Zachariasen [3-95].
        :return: PsiH.
        """
        return self._psi_H

    def PsiHBar(self):
        """
        Returns Psi0 as defined in Zachariasen [3-95].
        :return: PsiHBar.
        """
        return self._psi_H_bar

    def thickness(self):
        """
        Returns crystal thickness.
        :return: Thickness of the crystal.
        """
        return self._thickness

    def dSpacing(self):
        """
        Returns distance between the reflection planes.
        :return: Distance between the reflection planes.
        """
        return self._d_spacing

    def geometryType(self):
        """
        Returns the geometry types, i.e. BraggTransmission, LaueDiffraction,...
        :return: Geometry type.
        """
        return self._geometryType

    def log(self, string):
        """
        Logs a string.
        :param string: String to log.
        """
        print(string)

    def logDebug(self, string):
        """
        Logs a debug string.
        :param string: String to log.
        """
        self.log("<DEBUG>: " + string)

    def _calculateGamma(self, photon):
        """
        Calculates the projection cosine gamma as defined in Zachariasen [3-115].
        :param photon: Photon that is projected onto the surface normal.
        :return: Projection cosine gamma.
        """
        gamma = photon.unitDirectionVector().scalarProduct(self.surface_normal().getNormalizedVector())
        # Our crystal normal is pointing outside the crystal medium. Zachariasen's normal points
        # into the crystal medium (pag 112). Therefore, we change the sign.
        gamma = -gamma
        return gamma

    def _calculatePhotonOut(self, photon_in):
        """
        Solves the Laue equation to calculates the outgoing photon from the incoming photon and the Bragg normal.
        :param photon_in: Incoming photon or photon bunch.
        :return: Outgoing photon or photon bunch
        """
        # # Retrieve k_0.
        # k_in = photon_in.wavevector()

        # # Solve unscaled Laue equation.
        # k_out = self.braggNormal().addVector(k_in)

        # Create photon in k_out direction and scale by setting the photon energy.
        # photon_out = Photon(photon_in.energy(), k_out)
        """
        GENERAL VERSION:
        Solves the Laue equation for the parallel components of the vectors and
        uses the conservation of the wavevector modulus to calculate the outgoing wavevector
        even for diffraction not at the Bragg angle.
        """
        # Retrieve k_0.
        k_in = photon_in.wavevector()

        # Decompose the vector into a component parallel to the surface normal and
        # a component parallel to the surface: (k_in * n) n.
        k_in_normal = self.surface_normal().scalarMultiplication(k_in.scalarProduct(self.surface_normal()))
        k_in_parallel = k_in.subtractVector(k_in_normal)

        # Retrieve the B_H vector.
        B_H = self.braggNormal()

        # Decompose the vector into a component parallel to the surface normal and
        # a component parallel to the surface: (B_H * n) n.
        B_H_normal = self.surface_normal().scalarMultiplication(B_H.scalarProduct(self.surface_normal()))
        B_H_parallel = B_H.subtractVector(B_H_normal)

        # Apply the Laue formula for the parallel components.
        k_out_parallel = k_in_parallel.addVector(B_H_parallel)

        # Calculate K_out normal.
        k_out_normal_modulus = numpy.sqrt(k_in.norm() ** 2 - k_out_parallel.norm() ** 2)
        k_out_normal = self.surface_normal().scalarMultiplication(k_out_normal_modulus)

        # # Calculate the outgoing photon.
        # # changed srio@esrf.eu to negative normal component to take into account that crystal normal points
        # # outsize
        # k_out_1 = k_out_parallel.addVector(k_out_normal)
        # k_out_2 = k_out_parallel.scalarMultiplication(-1.0).addVector(k_out_normal)
        #
        # # select k_out_1 or k_out_2
        #
        # k_out_Ewald = photon_in.unitDirectionVector().scalarMultiplication(photon_in.wavenumber())
        # k_out_Ewald = k_out_Ewald.addVector(B_H)
        # k_out_Ewald = k_out_Ewald.getNormalizedVector()
        #
        # tmp1 = k_out_1.scalarProduct(k_out_Ewald)
        # tmp2 = k_out_2.scalarProduct(k_out_Ewald)

        # TODO: try to put some logic here
        if self.geometryType() == BraggDiffraction():
            k_out = k_out_parallel.addVector(k_out_normal)
        elif self.geometryType() == LaueDiffraction():
            k_out = k_out_parallel.addVector(k_out_normal.scalarMultiplication(-1.0))
        elif self.geometryType() == BraggTransmission():
            k_out = k_out_parallel.addVector(k_out_normal)
        elif self.geometryType() == LaueTransmission():
            k_out = k_out_parallel.addVector(k_out_normal.scalarMultiplication(-1.0))
        else:
            raise Exception


        photon_out = photon_in.duplicate()
        photon_out.setUnitDirectionVector(k_out)

        if self.isDebug:
            self.logDebug("surface normal" + str(self.surface_normal().components()))
            self.logDebug("Angle bragg normal photon_in"
                          + str((photon_in.unitDirectionVector().angle(self.braggNormal()),
                                numpy.pi * 0.5 - photon_in.unitDirectionVector().angle(self.braggNormal()))))
            self.logDebug("Angle bragg normal photon_out"
                          + str((photon_out.unitDirectionVector().angle(self.braggNormal()),
                                numpy.pi * 0.5 - photon_out.unitDirectionVector().angle(self.braggNormal()))))
            self.logDebug("photon_in direction" + str(photon_in.unitDirectionVector().components()))
            self.logDebug("photon_out direction" + str(photon_out.unitDirectionVector().components()))

        # Return outgoing photon.
        return photon_out


    def _calculateAlphaZac(self, photon_in):
        """
        Calculates alpha ("refraction index difference between waves in the crystal") as defined in Zachariasen [3-114b].
        :param photon_in: Incoming photon.
        :return: alpha.
        """
        # Calculate scalar product k_0 and B_H.
        k_0_times_B_h = photon_in.wavevector().scalarProduct(self.braggNormal())

        # Get norm k_0.
        wavenumber = photon_in.wavenumber()

        # Calculate alpha.
        zac_alpha = (wavenumber ** -2) * (self.braggNormal().norm() ** 2 + 2 * k_0_times_B_h)
        # we defined theta = theta_b + deviation (for symmetric Bragg), therefore the
        # approximated value of zac_alpha (compare eqs 3.116)  = 2 (theta_b-theta) sin(2 theta_b) =
        # = 2 (-deviation) * sin(2 theta_b)

        # Return alpha.
        return zac_alpha

    def _calculateGuigayAlpha(self, photon_in):
        """
        Calculates alpha ("refraction index difference between waves in the crystal") as defined in Zachariasen [3-114b].
        :param photon_in: Incoming photon.
        :return: alpha.
        """
        k0_dot_H = photon_in.wavevector().scalarProduct(self.braggNormal()) # scalar product k0 and H.
        wavenumber = photon_in.wavenumber() #  norm of k0.
        alpha = - (wavenumber ** -2) * (self.braggNormal().norm() ** 2 + 2 * k0_dot_H)
        return alpha


    def _calculateZacB(self, photon_in, photon_out):
        """
        Calculates asymmetry ratio b as defined in Zachariasen [3-115].
        :param photon_in: Incoming photon.
        :param photon_out: Outgoing photon.
        :return: Asymmetry ratio b.
        """
        # TODO: revise this algorithm, it is not exactly as in Zachariasen [3-115]
        numerator   = self.surface_normal().scalarProduct(photon_in.wavevector())
        denominator = self.surface_normal().scalarProduct(photon_out.wavevector())
        zac_b = numerator / denominator

        return zac_b


    def _calculateGuigayB(self, photon_in):
        """
        Calculates asymmetry ratio b as defined in Guigay.
        :param photon_in: Incoming photon.
        :return: Asymmetry ratio b.
        Note that this b changes when K0 (photon_in.wavevector()) changes
        """
        KH = photon_in.wavevector().addVector(self.braggNormal())
        if isinstance(photon_in, Photon):
            photon_outG = Photon(energy_in_ev=photon_in.energy(), direction_vector=KH)
            return self._calculateGamma(photon_in) / self._calculateGamma(photon_outG)
        elif isinstance(photon_in, PhotonBunch):
            photon_outG = PhotonBunch.initialize_from_energies_and_directions(photon_in.energy(), KH)
            return self._calculateGamma(photon_in) / self._calculateGamma(photon_outG)


    def _calculateZacQ(self, zac_b, effective_psi_h, effective_psi_h_bar):
        """
        Calculates q as defined in Zachariasen [3-123].
        :param zac_b: Asymmetry ratio b as defined in Zachariasen [3-115].
        :param effective_psi_h: Effective PsiH (depending of polarisation. See text following [3.-139]).
        :param effective_psi_h_bar: Effective PsiHBar (depending of polarisation. See text following [3.-139]).
        :return: q.
        """
        return zac_b * effective_psi_h * effective_psi_h_bar

    def _calculateZacZ(self, zac_b, zac_alpha):
        """
        Calcualtes z as defined in Zachariasen [3-123].
        :param zac_b: Asymmetry ratio b as defined in Zachariasen [3-115].
        :param zac_alpha: Diffraction index difference of crystal fields.
        :return: z.
        """
        return (1.0e0 - zac_b) * 0.5e0 * self.Psi0() + zac_b * 0.5e0 * zac_alpha

    def _createVariable(self, initial_value):
        """
        Factory method for calculation variable. Delegates to active calculation strategy.
        :param initial_value: Inital value of the variable.
        :return: Variable to use for the calculation.
        """
        return self._calculation_strategy.createVariable(initial_value)

    def _exponentiate(self, power):
        """
        Exponentiates to the power using active calculation strategy. (plain python or arbitrary precision)
        :param power: Calculation variable.
        :return: Exponential.
        """
        return self._calculation_strategy.exponentiate(self._createVariable(power))

    def _sin(self, power):
        """
        Sin to the power using active calculation strategy. (plain python or arbitrary precision)
        :param power: Calculation variable.
        :return: Sin.
        """
        return self._calculation_strategy.sin(self._createVariable(power))

    def _cos(self, power):
        """
        Cos to the power using active calculation strategy. (plain python or arbitrary precision)
        :param power: Calculation variable.
        :return: Cos.
        """
        return self._calculation_strategy.cos(self._createVariable(power))

    def _toComplex(self, variable):
        """
        Converts calculation variable to complex. Delegates to active calculation strategy.
        :param variable: Calculation variable.
        :return: Calculation variable as complex.
        """
        return self._calculation_strategy.toComplex(variable)

    def _calculateComplexAmplitude(self, photon_in, zac_q, zac_z, gamma_0, effective_psi_h_bar):
        """
        Calculates the complex amplitude of the questioned wave: diffracted or transmission.
        :param photon_in: Incoming photon.
        :param zac_q: q as defined in Zachariasen [3-123].
        :param zac_z: z as defined in Zachariasen [3-123].
        :param gamma_0: Projection cosine as defined in Zachariasen [3-115].
        :param effective_psi_h_bar: Effective PsiHBar (depending of polarisation. See text following [3.-139]).
        :return: Complex amplitude.
        """
        # Calculate geometry independent parts.
        tmp_root = (zac_q + zac_z * zac_z) ** 0.5

        zac_x1 = (-1.0 * zac_z + tmp_root) / effective_psi_h_bar
        zac_x2 = (-1.0 * zac_z - tmp_root) / effective_psi_h_bar
        zac_delta1 = 0.5 * (self.Psi0() - zac_z + tmp_root)
        zac_delta2 = 0.5 * (self.Psi0() - zac_z - tmp_root)
        zac_phi1 = 2 * numpy.pi / gamma_0 / photon_in.wavelength() * zac_delta1
        zac_phi2 = 2 * numpy.pi / gamma_0 / photon_in.wavelength() * zac_delta2
       
        zac_c1 = -1j * self.thickness() * zac_phi1
        zac_c2 = -1j * self.thickness() * zac_phi2

        if self.isDebug:
            self.logDebug("__zac_c1" + str(zac_c1))
            self.logDebug("__zac_c2" + str(zac_c2))

        cv_zac_c1 = self._exponentiate(zac_c1)
        cv_zac_c2 = self._exponentiate(zac_c2)

        cv_zac_x1 = self._createVariable(zac_x1)
        cv_zac_x2 = self._createVariable(zac_x2)

        # Calculate complex amplitude according to given geometry.
        if self.geometryType() == BraggDiffraction():
            complex_amplitude = cv_zac_x1 * cv_zac_x2 * (cv_zac_c2 - cv_zac_c1) / \
                                (cv_zac_c2 * cv_zac_x2 - cv_zac_c1 * cv_zac_x1)
        elif self.geometryType() == LaueDiffraction():
            complex_amplitude = cv_zac_x1 * cv_zac_x2 * (cv_zac_c1 - cv_zac_c2) / \
                                (cv_zac_x2 - cv_zac_x1)
        elif self.geometryType() == BraggTransmission():
            complex_amplitude = cv_zac_c1 * cv_zac_c2 * (cv_zac_x2 - cv_zac_x1) / \
                                (cv_zac_c2 * cv_zac_x2 - cv_zac_c1 * cv_zac_x1)
        elif self.geometryType() == LaueTransmission():
            complex_amplitude = (cv_zac_x2 * cv_zac_c1 - cv_zac_x1 * cv_zac_c2) / \
                                (cv_zac_x2 - cv_zac_x1)
        else:
            raise Exception

        if self.isDebug:
            self.logDebug("ctemp: " + str(tmp_root))
            self.logDebug("zac_z" + str(zac_z))
            self.logDebug("zac_q" + str(zac_q))
            self.logDebug("zac delta 1" + str(zac_delta1))
            self.logDebug("zac delta 2" + str(zac_delta2))
            self.logDebug("gamma_0" + str(gamma_0))
            self.logDebug("wavelength" + str(photon_in.wavelength()))
            self.logDebug("zac phi 1" + str(zac_phi1))
            self.logDebug("zac phi 2" + str(zac_phi2))
            self.logDebug("zac_c1: " + str(cv_zac_c1))
            self.logDebug("zac_c2: " + str(cv_zac_c2))
            self.logDebug("zac_x1: " + str(cv_zac_x1))
            self.logDebug("zac_x2: " + str(cv_zac_x2))

        # return ComplexAmplitude(complex(complex_amplitude))
        return complex_amplitude # ComplexAmplitude(complex_amplitude)

    def _calculatePolarizationS(self, photon_in, zac_b, zac_z, gamma_0):
        """
        Calculates complex amplitude for the S polarization.
        :param photon_in: Incoming photon.
        :param zac_z: z as defined in Zachariasen [3-123].
        :param gamma_0: Projection cosine as defined in Zachariasen [3-115].
        :return: Complex amplitude of S polarization.
        """
        zac_q = self._calculateZacQ(zac_b, self.PsiH(), self.PsiHBar())

        return self._calculateComplexAmplitude(photon_in, zac_q, zac_z, gamma_0,
                                               self.PsiHBar())

    def _calculatePolarizationP(self, photon_in, zac_b, zac_z, gamma_0):
        """
        Calculates complex amplitude for the P polarization.
        :param photon_in: Incoming photon.
        :param zac_b: Asymmetry ratio b as defined in Zachariasen [3-115].
        :param zac_z: z as defined in Zachariasen [3-123].
        :param gamma_0: Projection cosine as defined in Zachariasen [3-115].
        :return: Complex amplitude of P polarization.
        """
        effective_psi_h = self.PsiH() * numpy.cos(2 * self.braggAngle())
        effective_psi_h_bar = self.PsiHBar() * numpy.cos(2 * self.braggAngle())

        zac_q = self._calculateZacQ(zac_b, effective_psi_h, effective_psi_h_bar)

        return self._calculateComplexAmplitude(photon_in, zac_q, zac_z, gamma_0,
                                               effective_psi_h_bar)

    def calculateDiffraction(self,
                             photon_in,
                             calculation_method=0, # 0=Zachariasen, 1=Guigay
                             is_thick=0, # for Guigay only
                             use_transfer_matrix=0, # for Guigay only
                             ):
        """
        Calculate diffraction for incoming photon.
        :param photon_in: Incoming photon or Photon bunch.
        :return: Complex amplitude of the diffraction.
        """

        if calculation_method == 0:
            # print(">>>> Using Zachariasen equations...")
            return self.calculateDiffractionZachariasen(photon_in)
        else:
            # print(">>>> Using Guigay equations...")
            return self.calculateDiffractionGuigay(photon_in, is_thick=is_thick, use_transfer_matrix=use_transfer_matrix)


    def calculateDiffractionZachariasen(self, photon_in):
        """
        Calculate diffraction for incoming photon.
        :param photon_in: Incoming photon.
        :return: Complex amplitude of the diffraction.
        """
        # Initialize return variable.

        result = {"S": None,
                  "P": None}


        # Calculate photon out.
        photon_out = self._calculatePhotonOut(photon_in)

        # Calculate crystal field refraction index difference.
        zac_alpha = self._calculateAlphaZac(photon_in)

        # Calculate asymmetry ratio.
        # zac_b = self._calculateZacB(photon_in, photon_out)
        zac_b = self._calculateGuigayB(photon_in)  # todo: check if this is the same for Zac

        # Calculate z as defined in Zachariasen [3-123].
        zac_z = self._calculateZacZ(zac_b, zac_alpha)

        # Calculate projection cosine.
        gamma_0 = self._calculateGamma(photon_in)

        # Calculate complex amplitude for S and P polarization.
        result["S"] = self._calculatePolarizationS(photon_in, zac_b, zac_z, gamma_0)
        result["P"] = self._calculatePolarizationP(photon_in, zac_b, zac_z, gamma_0)

        # Note division by |b| in intensity (thus sqrt(|b|) in amplitude)
        # for power balance (see Zachariasen pag. 122)
        #
        # This factor only applies to diffracted beam, not to transmitted beams
        # (see private communication M. Rio (ESRF) and J. Sutter (DLS))
        if (self.geometryType() == BraggDiffraction() or
                self.geometryType() == LaueDiffraction()):
            # result["S"].rescale(1.0 / sqrt(abs(zac_b)))
            # result["P"].rescale(1.0 / sqrt(abs(zac_b)))
            result["S"] *= (1.0 / numpy.sqrt(abs(zac_b)))
            result["P"] *= (1.0 / numpy.sqrt(abs(zac_b)))

        # If debugging output is turned on.
        if self.isDebug:
            self._logMembers(zac_b, zac_alpha, photon_in, photon_out, result)

        # Returns the complex amplitudes.
        return result

    def calculateDiffractionGuigay(self, photon_in, debug=0, s_ratio=None,
                                   is_thick=0,
                                   use_transfer_matrix=0, # is faster to use use_transfer_matrix=0
                                   ):
        """
        Calculate diffraction for incoming photon.
        :param photon_in: Incoming photon.
        :return: Complex amplitude of the diffraction.
        """
        # Initialize return variable.

        result = {"S": None,
                  "P": None}

        # if isinstance(photon_in, Photon):
        #     result = {}
        # elif isinstance(photon_in, PhotonBunch):
        #     result = [{}] * photon_in.getNumberOfPhotons()

        if use_transfer_matrix:
            transfer_matrix_s = self.calculateTransferMatrix(photon_in, polarization=0, is_thick=is_thick)
            m11_s, m12_s, m21_s, m22_s = transfer_matrix_s
            scattering_matrix_s = self.calculateScatteringMatrixFromTransferMatrix(transfer_matrix_s)
            s11_s, s12_s, s21_s, s22_s = scattering_matrix_s

            transfer_matrix_p = self.calculateTransferMatrix(photon_in, polarization=1, is_thick=is_thick)
            m11_p, m12_p, m21_p, m22_p = transfer_matrix_p
            scattering_matrix_p = self.calculateScatteringMatrixFromTransferMatrix(transfer_matrix_p)
            s11_p, s12_p, s21_p, s22_p = scattering_matrix_p

            # result["transfer_matrix_s"] = transfer_matrix_s
            # result["transfer_matrix_p"] = transfer_matrix_p
            # result["scattering_matrix_s"] = scattering_matrix_s
            # result["scattering_matrix_p"] = scattering_matrix_p
            result["m11_s"] = m11_s
            result["m12_s"] = m12_s
            result["m21_s"] = m21_s
            result["m22_s"] = m22_s
            result["m11_p"] = m11_p
            result["m12_p"] = m12_p
            result["m21_p"] = m21_p
            result["m22_p"] = m22_p

            result["s11_s"] = s11_s
            result["s12_s"] = s12_s
            result["s21_s"] = s21_s
            result["s22_s"] = s22_s
            result["s11_p"] = s11_p
            result["s12_p"] = s12_p
            result["s21_p"] = s21_p
            result["s22_p"] = s22_p

            if self.geometryType() == BraggDiffraction():
                # guigay, sanchez del rio,  eq 31a
                complex_amplitude_s = s21_s
                complex_amplitude_p = s21_p
            elif self.geometryType() == BraggTransmission():
                # guigay, sanchez del rio,  eq 31b
                complex_amplitude_s = s11_s
                complex_amplitude_p = s11_p
            elif self.geometryType() == LaueDiffraction():
                # guigay, sanchez del rio,  eq 27a
                complex_amplitude_s = m21_s
                complex_amplitude_p = m21_p
            elif self.geometryType() == LaueTransmission():
                # guigay, sanchez del rio,  eq 27b
                complex_amplitude_s = m11_s
                complex_amplitude_p = m11_p
            else:
                raise Exception

        else:  # todo: try to fix this and then delete........
            # Calculate photon out.
            photon_out = self._calculatePhotonOut(photon_in)

            # Calculate crystal field refraction index difference.
            # Note that Guigay's definition of alpha has the opposite sign as in Zachariasen!
            alpha = self._calculateGuigayAlpha(photon_in)
            if debug: print("guigay alpha: ", alpha)

            guigay_b = self._calculateGuigayB(photon_in) # gamma_0 / gamma_H
            if debug: print("guigay_b: ", guigay_b)

            gamma_0 = self._calculateGamma(photon_in)
            T = self.thickness() / gamma_0
            if debug:
                print("gamma_0: ", gamma_0)
                print("T: ", T)

            effective_psi_0 = numpy.conjugate(self.Psi0())  # I(Psi0) > 0 (for absorption!!)

            w = guigay_b * (alpha / 2) + effective_psi_0 * (guigay_b - 1) / 2
            omega = numpy.pi / photon_in.wavelength() * w

            if self.geometryType() == BraggDiffraction():
                if s_ratio is None:
                    s = 0.0
                else:
                    s = T * s_ratio
                # sigma polarization
                effective_psi_h = numpy.conjugate(self.PsiH())
                effective_psi_h_bar = numpy.conjugate(self.PsiHBar())
                uh = effective_psi_h * numpy.pi / photon_in.wavelength()
                uh_bar = effective_psi_h_bar * numpy.pi / photon_in.wavelength()
                u0 = effective_psi_0 * numpy.pi / photon_in.wavelength()

                # guigay, sanchez del rio,  eq 31a
                if is_thick == 0:
                    SQ = numpy.sqrt(guigay_b * effective_psi_h * effective_psi_h_bar + w ** 2)
                    a = numpy.pi / photon_in.wavelength() * SQ
                    complex_amplitude_s = 1j * guigay_b * uh * self._sin(a * s - a * T) / \
                                        (a * self._cos(a * T) + 1j * omega * self._sin(a * T)) * \
                                        self._exponentiate(1j * s * (omega + u0))
                else:
                    #Thickkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkk!
                    asquared = (numpy.pi / photon_in.wavelength())**2 * (guigay_b * effective_psi_h * effective_psi_h_bar + w ** 2)
                    aa = 1 / numpy.sqrt(2) * ( (asquared).imag / numpy.sqrt(numpy.abs(asquared)-(asquared).real) + \
                                               1j * numpy.sqrt(numpy.abs(asquared) - (asquared).real))

                    complex_amplitude_s = (aa + omega) / uh_bar

                # print(">>> complex_amplitude_s: ", complex_amplitude_s)
                # pi polarization
                effective_psi_h = numpy.conjugate(self.PsiH()) * numpy.cos(2 * self.braggAngle())
                effective_psi_h_bar = numpy.conjugate(self.PsiHBar()) * numpy.cos(2 * self.braggAngle())
                uh = effective_psi_h * numpy.pi / photon_in.wavelength()
                uh_bar = effective_psi_h_bar * numpy.pi / photon_in.wavelength()
                u0 = effective_psi_0 * numpy.pi / photon_in.wavelength()

                # guigay, sanchez del rio,  eq 31b
                if is_thick == 0:
                    SQ = numpy.sqrt(guigay_b * effective_psi_h * effective_psi_h_bar + w ** 2)
                    a = numpy.pi / photon_in.wavelength() * SQ
                    complex_amplitude_p = 1j * guigay_b * uh * self._sin( a * s - a * T) / \
                                        (a * self._cos(a * T) + 1j * omega * self._sin(a * T)) * \
                                        self._exponentiate(1j * s * (omega + u0))
                else:
                    #Thickkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkk!
                    asquared = (numpy.pi / photon_in.wavelength())**2 * (guigay_b * effective_psi_h * effective_psi_h_bar + w ** 2)
                    aa = 1 / numpy.sqrt(2) * ( (asquared).imag / numpy.sqrt(numpy.abs(asquared)-(asquared).real) + \
                                               1j * numpy.sqrt(numpy.abs(asquared) - (asquared).real))

                    complex_amplitude_p = (aa + omega) / uh_bar

            elif self.geometryType() == BraggTransmission():
                if s_ratio is None:
                    s = T
                else:
                    s = T * s_ratio
                # sigma polarization
                effective_psi_h = numpy.conjugate(self.PsiH())
                effective_psi_h_bar = numpy.conjugate(self.PsiHBar())
                uh_bar = effective_psi_h_bar * numpy.pi / photon_in.wavelength()
                u0 = effective_psi_0 * numpy.pi / photon_in.wavelength()


                # guigay, sanchez del rio,  eq 31b
                # complex_amplitude_s = 1 / (a * self._cos(a * T) + 1j * omega * self._sin(a * T))
                # complex_amplitude_s *= numpy.exp(1j * T * (omega + u0))
                if is_thick == 0:
                    SQ = numpy.sqrt(guigay_b * effective_psi_h * effective_psi_h_bar + w ** 2)
                    a = numpy.pi / photon_in.wavelength() * SQ

                    complex_amplitude_s = (a * self._cos(a * s - a * T) - 1j * omega * self._sin(a * s - a * T)) / \
                                          (a * self._cos(a * T) + 1j * omega * self._sin(a * T))
                    complex_amplitude_s *= numpy.exp(1j * T * (omega + u0))
                else:
                    #Thickkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkk!
                    asquared = (numpy.pi / photon_in.wavelength())**2 * (guigay_b * effective_psi_h * effective_psi_h_bar + w ** 2)
                    aa = 1 / numpy.sqrt(2) * ( (asquared).imag / numpy.sqrt(numpy.abs(asquared)-(asquared).real) + \
                                               1j * numpy.sqrt(numpy.abs(asquared) - (asquared).real))

                    complex_amplitude_s = 2 * aa / (aa - omega) * numpy.exp(1j * T * (u0 + omega + aa))

                # pi polarization
                effective_psi_h = numpy.conjugate(self.PsiH()) * numpy.cos(2 * self.braggAngle())
                effective_psi_h_bar = numpy.conjugate(self.PsiHBar()) * numpy.cos(2 * self.braggAngle())
                u0 = effective_psi_0 * numpy.pi / photon_in.wavelength()

                # guigay, sanchez del rio,  eq 31b
                # complex_amplitude_p = 1 / (a * self._cos(a * T) + 1j * omega * self._sin(a * T))
                # complex_amplitude_p *= numpy.exp(1j * T * (omega + u0))
                if is_thick == 0:
                    SQ = numpy.sqrt(guigay_b * effective_psi_h * effective_psi_h_bar + w ** 2)
                    a = numpy.pi / photon_in.wavelength() * SQ

                    complex_amplitude_p = (a * self._cos(a * s - a * T) - 1j * omega * self._sin(a * s - a * T)) / \
                                          (a * self._cos(a * T) + 1j * omega * self._sin(a * T))
                    complex_amplitude_p *= numpy.exp(1j * T * (omega + u0))
                else:
                    # Thickkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkk!
                    asquared = (numpy.pi / photon_in.wavelength()) ** 2 * (
                                guigay_b * effective_psi_h * effective_psi_h_bar + w ** 2)
                    aa = 1 / numpy.sqrt(2) * ((asquared).imag / numpy.sqrt(numpy.abs(asquared) - (asquared).real) + \
                                              1j * numpy.sqrt(numpy.abs(asquared) - (asquared).real))

                    complex_amplitude_p = 2 * aa / (aa - omega) * numpy.exp(1j * T * (u0 + omega + aa))

            elif self.geometryType() == LaueDiffraction():
                if s_ratio is None:
                    s = T
                else:
                    s = T * s_ratio

                # sigma polarization
                effective_psi_h     = numpy.conjugate(self.PsiH())
                effective_psi_h_bar = numpy.conjugate(self.PsiHBar())
                uh = effective_psi_h * numpy.pi / photon_in.wavelength()
                u0 = effective_psi_0 * numpy.pi / photon_in.wavelength()

                # guigay, sanchez del rio,  eq 27a todo: as a function of s
                if is_thick == 0:
                    SQ = numpy.sqrt(guigay_b * effective_psi_h * effective_psi_h_bar + w ** 2)
                    a = numpy.pi / photon_in.wavelength() * SQ
                    complex_amplitude_s = 1j * guigay_b * uh * self._sin(a * s) / a
                    complex_amplitude_s *= self._exponentiate(1j * s * (omega + u0))
                else:
                    # Thickkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkk!
                    asquared = (numpy.pi / photon_in.wavelength()) ** 2 * (
                                guigay_b * effective_psi_h * effective_psi_h_bar + w ** 2)
                    aa = 1 / numpy.sqrt(2) * ((asquared).imag / numpy.sqrt(numpy.abs(asquared) - (asquared).real) + \
                                              1j * numpy.sqrt(numpy.abs(asquared) - (asquared).real))

                    complex_amplitude_s = - guigay_b * uh / (2 * aa) * self._exponentiate(1j * s * (omega + u0 - aa))


                # pi polarization
                effective_psi_h     = numpy.conjugate(self.PsiH()) * numpy.cos(2 * self.braggAngle())
                effective_psi_h_bar = numpy.conjugate(self.PsiHBar()) * numpy.cos(2 * self.braggAngle())
                uh = effective_psi_h * numpy.pi / photon_in.wavelength()
                u0 = effective_psi_0 * numpy.pi / photon_in.wavelength()

                # guigay, sanchez del rio,  eq 27a todo: as a function of s
                if is_thick == 0:
                    SQ = numpy.sqrt(guigay_b * effective_psi_h * effective_psi_h_bar + w ** 2)
                    a = numpy.pi / photon_in.wavelength() * SQ
                    complex_amplitude_p = 1j * guigay_b * uh * self._sin(a * s) / a
                    complex_amplitude_p *= self._exponentiate(1j * s * (omega + u0))
                else:
                    # Thickkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkk!
                    asquared = (numpy.pi / photon_in.wavelength()) ** 2 * (
                                guigay_b * effective_psi_h * effective_psi_h_bar + w ** 2)
                    aa = 1 / numpy.sqrt(2) * ((asquared).imag / numpy.sqrt(numpy.abs(asquared) - (asquared).real) + \
                                              1j * numpy.sqrt(numpy.abs(asquared) - (asquared).real))

                    complex_amplitude_p = - guigay_b * uh / (2 * aa) * self._exponentiate(1j * s * (omega + u0 - aa))


            elif self.geometryType() == LaueTransmission():
                if s_ratio is None:
                    s = T
                else:
                    s = T * s_ratio

                # sigma polarization
                effective_psi_h = numpy.conjugate(self.PsiH())
                effective_psi_h_bar = numpy.conjugate(self.PsiHBar())
                u0 = effective_psi_0 * numpy.pi / photon_in.wavelength()


                # guigay, sanchez del rio,  eq 27b todo: as a function of s
                if is_thick == 0:
                    SQ = numpy.sqrt(guigay_b * effective_psi_h * effective_psi_h_bar + w ** 2)
                    a = numpy.pi / photon_in.wavelength() * SQ
                    complex_amplitude_s = numpy.cos(a * s) - 1j * omega * self._sin(a * s) / a
                    complex_amplitude_s *= self._exponentiate(1j * s * (omega + u0))
                else:
                    # Thickkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkk!
                    asquared = (numpy.pi / photon_in.wavelength()) ** 2 * (
                                guigay_b * effective_psi_h * effective_psi_h_bar + w ** 2)
                    aa = 1 / numpy.sqrt(2) * ((asquared).imag / numpy.sqrt(numpy.abs(asquared) - (asquared).real) + \
                                              1j * numpy.sqrt(numpy.abs(asquared) - (asquared).real))

                    complex_amplitude_s = self._exponentiate(1j * s * (omega + u0 - aa)) * 0.5 * (1 + omega / aa)

                # pi polarization
                effective_psi_h = numpy.conjugate(self.PsiH()) * numpy.cos(2 * self.braggAngle())
                effective_psi_h_bar = numpy.conjugate(self.PsiHBar()) * numpy.cos(2 * self.braggAngle())
                u0 = effective_psi_0 * numpy.pi / photon_in.wavelength()

                # guigay, sanchez del rio,  eq 27b todo: as a function of s
                if is_thick == 0:
                    SQ = numpy.sqrt(guigay_b * effective_psi_h * effective_psi_h_bar + w ** 2)
                    a = numpy.pi / photon_in.wavelength() * SQ
                    complex_amplitude_p = numpy.cos(a * s) - 1j * omega * self._sin(a * s) / a
                    complex_amplitude_p *= self._exponentiate(1j * s * (omega + u0))
                else:
                    # Thickkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkk!
                    asquared = (numpy.pi / photon_in.wavelength()) ** 2 * (
                                guigay_b * effective_psi_h * effective_psi_h_bar + w ** 2)
                    aa = 1 / numpy.sqrt(2) * ((asquared).imag / numpy.sqrt(numpy.abs(asquared) - (asquared).real) + \
                                              1j * numpy.sqrt(numpy.abs(asquared) - (asquared).real))

                    complex_amplitude_p = self._exponentiate(1j * s * (omega + u0 - aa)) * 0.5 * (1 + omega / aa)

            else:
                raise Exception

            result["alpha"] = alpha
            result["b"] = guigay_b

        # Calculate complex amplitude for S and P polarization.


        if isinstance(photon_in, Photon):
            result["S"] = complex_amplitude_s # ComplexAmplitude(complex(complex_amplitude_s))
            result["P"] = complex_amplitude_p # ComplexAmplitude(complex(complex_amplitude_p))
            result["s"] = complex_amplitude_s # these coefficients will not be weighted for power.
            result["p"] = complex_amplitude_p # these coefficients will not be weighted for power.
            # Note division by |b| in intensity (thus sqrt(|b|) in amplitude)
            # for power balance (see Zachariasen pag. 122)
            #
            # This factor only applies to diffracted beam, not to transmitted beams
            # (see private communication M. Rio (ESRF) and J. Sutter (DLS))

            if (self.geometryType() == BraggDiffraction() or self.geometryType() == LaueDiffraction()):
                # recalculate guigay_b ...
                photon_out = self._calculatePhotonOut(photon_in)
                # alpha = self._calculateGuigayAlpha(photon_in)
                guigay_b = self._calculateGuigayB(photon_in)  # gamma_0 / gamma_H
                # result["S"].rescale(1.0 / sqrt(abs(guigay_b)))
                # result["P"].rescale(1.0 / sqrt(abs(guigay_b)))
                result["S"] *= 1.0 / numpy.sqrt(abs(guigay_b))
                result["P"] *= 1.0 / numpy.sqrt(abs(guigay_b))

            # If debugging output is turned on.
            if self.isDebug:
                self._logMembers(guigay_b, zac_alpha, photon_in, photon_out, result)
        elif isinstance(photon_in, PhotonBunch):
            result["S"] = complex_amplitude_s # ComplexAmplitude(complex(complex_amplitude_s))
            result["P"] = complex_amplitude_p # ComplexAmplitude(complex(complex_amplitude_p))
            result["s"] = complex_amplitude_s
            result["p"] = complex_amplitude_p
            if (self.geometryType() == BraggDiffraction() or self.geometryType() == LaueDiffraction()):
                # recalculate guigay_b ...
                photon_out = self._calculatePhotonOut(photon_in)
                # alpha = self._calculateGuigayAlpha(photon_in)
                guigay_b = self._calculateGuigayB(photon_in)  # gamma_0 / gamma_H

                result["S"] *= 1.0 / numpy.sqrt(abs(guigay_b))
                result["P"] *= 1.0 / numpy.sqrt(abs(guigay_b))

            # for i in range(photon_in.getNumberOfPhotons()):
            #     result[i]["S"] = ComplexAmplitude(complex(complex_amplitude_s[i]))
            #     result[i]["P"] = ComplexAmplitude(complex(complex_amplitude_p[i]))
            #     result[i]["s"] = complex_amplitude_s[i]
            #     result[i]["p"] = complex_amplitude_p[i]
            # if (self.geometryType() == BraggDiffraction() or self.geometryType() == LaueDiffraction()):
            #     # recalculate guigay_b ...
            #     photon_out = self._calculatePhotonOut(photon_in)
            #     # alpha = self._calculateGuigayAlpha(photon_in)
            #     guigay_b = self._calculateGuigayB(photon_in)  # gamma_0 / gamma_H
            #     for i in range(photon_in.getNumberOfPhotons()):
            #         result[i]["S"].rescale(1.0 / sqrt(abs(guigay_b[i])))
            #         result[i]["P"].rescale(1.0 / sqrt(abs(guigay_b[i])))


        # # Note division by |b| in intensity (thus sqrt(|b|) in amplitude)
        # # for power balance (see Zachariasen pag. 122)
        # #
        # # This factor only applies to diffracted beam, not to transmitted beams
        # # (see private communication M. Rio (ESRF) and J. Sutter (DLS))
        #
        # if (self.geometryType() == BraggDiffraction() or self.geometryType() == LaueDiffraction()):
        #
        #     # recalculate guigay_b ...
        #     photon_out = self._calculatePhotonOut(photon_in)
        #     alpha = self._calculateGuigayAlpha(photon_in)
        #     guigay_b = self._calculateGuigayB(photon_in) # gamma_0 / gamma_H
        #
        #
        #     result["S"].rescale(1.0 / sqrt(abs(guigay_b)))
        #     result["P"].rescale(1.0 / sqrt(abs(guigay_b)))
        #
        # # If debugging output is turned on.
        # if self.isDebug:
        #     self._logMembers(guigay_b, zac_alpha, photon_in, photon_out, result)

        # Returns the complex amplitudes.
        return result

    def calculateTransferMatrix(self, photon_in, polarization=0, is_thick=0):

        photon_out = self._calculatePhotonOut(photon_in)
        alpha = self._calculateGuigayAlpha(photon_in)
        guigay_b = self._calculateGuigayB(photon_in)
        gamma_0 = self._calculateGamma(photon_in)
        T = self.thickness() / gamma_0

        if polarization == 0:
            pol_factor = 1.0
        else:
            pol_factor = numpy.cos(2 * self.braggAngle())

        effective_psi_0 = numpy.conjugate(self.Psi0())  # I(Psi0) > 0 (for absorption!!)

        w = guigay_b * (alpha / 2) + effective_psi_0 * (guigay_b - 1) / 2
        omega = numpy.pi / photon_in.wavelength() * w

        effective_psi_h     = numpy.conjugate(self.PsiH()) * pol_factor
        effective_psi_h_bar = numpy.conjugate(self.PsiHBar()) * pol_factor

        uh = effective_psi_h * numpy.pi / photon_in.wavelength()
        uh_bar = effective_psi_h_bar * numpy.pi / photon_in.wavelength()
        u0 = effective_psi_0 * numpy.pi / photon_in.wavelength()

        # guigay, sanchez del rio,  eq 26
        # phase_term = self._exponentiate(1j * T * (omega + u0))
        # m11 = self._cos(a * T) - 1j * omega * self._sin(a * T) / a
        # m12 = 1j *  uh_bar * self._sin(a * T) / a
        # m21 = 1j * guigay_b * uh * self._sin(a * T) / a
        # m22 = self._cos(a * T) + 1j * omega * self._sin(a * T) / a

        if is_thick:
            asquared = (numpy.pi / photon_in.wavelength()) ** 2 * (
                        guigay_b * effective_psi_h * effective_psi_h_bar + w ** 2)
            aa = 1 / numpy.sqrt(2) * ((asquared).imag / numpy.sqrt(numpy.abs(asquared) - (asquared).real) + \
                                      1j * numpy.sqrt(numpy.abs(asquared) - (asquared).real))

            phase_term = numpy.exp(1j * T * (omega + u0))
            sin_aT = 1j / 2 * self._exponentiate(-1j * aa * T) # self._sin(a * T)
            cos_aT = 1  / 2 * self._exponentiate(-1j * aa * T) # self._cos(a * T)
            # develop in series to avoid using mpmath (not working for Laue!!!!!!!!!!!!!!!)
            # x = aa * T
            # sin_aT = 1j / 2 * (1 -1j * x ) # self._sin(a * T)
            # cos_aT = 1  / 2 * (1 -1j * x ) # self._cos(a * T)


            m11 = cos_aT - 1j * omega * sin_aT / aa
            m12 = 1j *  uh_bar * sin_aT / aa
            m21 = 1j * guigay_b * uh * sin_aT / aa
            m22 = cos_aT + 1j * omega * sin_aT / aa
        else:
            SQ = numpy.sqrt(guigay_b * effective_psi_h * effective_psi_h_bar + w ** 2)
            a = numpy.pi / photon_in.wavelength() * SQ

            phase_term = numpy.exp(1j * T * (omega + u0))
            sin_aT = self._sin(a * T)
            cos_aT = self._cos(a * T)
            m11 = cos_aT - 1j * omega * sin_aT / a
            m12 = 1j *  uh_bar * sin_aT / a
            m21 = 1j * guigay_b * uh * sin_aT / a
            m22 = cos_aT + 1j * omega * sin_aT / a

        return m11 * phase_term, m12 * phase_term, m21 * phase_term, m22 * phase_term

    @classmethod
    def calculateScatteringMatrixFromTransferMatrix(self, transfer_matrix):
        m11, m12, m21, m22 = transfer_matrix
        s11 = m11 - m12 * m21 / m22
        s12 = m12 / m22
        s21 = -m21 / m22
        s22 = 1 / m22
        return s11, s12, s21, s22

    def calculateScatteringMatrix(self, photon_in, polarization=0):
        transfer_matrix = self.calculateTransferMatrix(photon_in, polarization=polarization)
        return self.calculateScatteringMatrixFromTransferMatrix(transfer_matrix)


    def _logMembers(self, zac_b, zac_alpha, photon_in, photon_out, result):
        """
        Debug logs the member variables and other relevant partial results.
        :param zac_b: Asymmetry ratio b
        :param zac_alpha: Diffraction index difference of crystal fields.
        :param photon_in: Incoming photon.
        :param result: Resulting complex amplitudes of the diffraction/transmission.
        """
        self.logDebug("Bragg angle: %f degrees \n" % (self.braggAngle() * 180 / pi))
        self.logDebug("psi0: (%.14f , %.14f)" % (self.Psi0().real, self.Psi0().imag))
        self.logDebug("psiH: (%.14f , %.14f)" % (self.PsiH().real, self.PsiH().imag))
        self.logDebug("psiHbar: (%.14f , %.14f)" % (self.PsiHBar().real, self.PsiHBar().imag))
        self.logDebug("d_spacing: %g " % self.dSpacing())
        self.logDebug('BraggNormal: ' + str(self.braggNormal().components()))
        self.logDebug('BraggNormal(Normalized): ' + str(self.braggNormal().getNormalizedVector().components()))
        self.logDebug('b(exact): ' + str(zac_b))
        self.logDebug('alpha: ' + str(zac_alpha))
        self.logDebug('k_0 wavelength: ' + str(photon_in.wavelength()))
        self.logDebug('PhotonInDirection:  ' + str(photon_in.unitDirectionVector().components()))
        self.logDebug('PhotonOutDirection: ' + str(photon_out.unitDirectionVector().components()))
        self.logDebug('comp ampl S: ' + str(result["S"].intensity()) + str(result["S"].phase()))
        self.logDebug('comp ampl P: ' + str(result["P"].intensity()) + str(result["P"].phase()))


if __name__ == "__main__":
    a = CalculationStrategyMPMath()

    a0 = a.createVariable(numpy.array(0))
    print("cos(0)", a.cos(a0))

    api = a.createVariable(numpy.array(numpy.pi))
    print("cos(pi)", a.cos(api))

    api = a.createVariable(numpy.array([numpy.pi] * 10))
    print("cos(pi)", a.cos(api))
    # print("exp(pi)", a.exponentiate(api))