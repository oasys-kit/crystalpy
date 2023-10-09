"""
Calculates crystal diffraction according to Guigay and Zachariasen formalisms of the dynamic theory of crystal diffraction
for perfect crystals.
Except for energy all units are in SI. Energy is in eV.
"""

import numpy
from crystalpy.util.Photon import Photon
from crystalpy.util.ComplexAmplitudePhotonBunch import ComplexAmplitudePhotonBunch
from crystalpy.diffraction.GeometryType import BraggDiffraction, LaueDiffraction, BraggTransmission, LaueTransmission

# Use mpmath if possible. Otherwise use numpy.
try:
    # raise ImportError
    import mpmath
    use_mpmath = True
    mpmath_sin = numpy.vectorize(mpmath.sin)
    mpmath_cos = numpy.vectorize(mpmath.cos)
    mpmath_exp = numpy.vectorize(mpmath.exp)
except ImportError:
    use_mpmath = False
    print("mpmath module for arbitrary-precision floating-point arithmetic could not be found!\n "
          "Using numpy instead. This could lead to overflow errors.\n")

if use_mpmath:
    print("Using mpmath.")
else:
    print("Using numpy instead of mpmath (not available).")


class CalculationStrategy(object):
    """Abstract strategy for calculation. Can be plain python or arbitrary precision like mpmath."""
    def createVariable(self, initial_value):
        """Factory method for calculation variable.

        Parameters
        ----------
        initial_value :
            Initial value of the variable.

        Returns
        -------
        type
            Calculation variable.

        """
        raise Exception("Must override this method.")

    def exponentiate(self, power):
        """Exponentiates to the power.

        Parameters
        ----------
        power :
            The power to raise to.

        Returns
        -------
        type
            Exponential.

        """
        raise Exception("Must override this method.")

    def sin(self, power):
        """Sin to the power.

        Parameters
        ----------
        power :
            The power to raise to.

        Returns
        -------
        type
            Sin.

        """
        raise Exception("Must override this method.")

    def cos(self, power):
        """Cos to the power.

        Parameters
        ----------
        power :
            The power to raise to.

        Returns
        -------
        type
            Cos.

        """
        raise Exception("Must override this method.")


    def toComplex(self, variable):
        """Converts calculation variable to native python complex.

        Parameters
        ----------
        variable :
            Calculation variable to convert.

        Returns
        -------
        type
            Native python complex variable.

        """
        raise Exception("Must override this method.")


class CalculationStrategyMPMath(CalculationStrategy):
    """Use mpmath for calculation."""
    def __init__(self):
        """
        Constructor.
        """
        # Use 32 digits in mpmath calculations.
        mpmath.mp.dps = 32

    def createVariable(self, initial_value):
        """Factory method for calculation variable.

        Parameters
        ----------
        initial_value :
            Initial value of the variable.

        Returns
        -------
        type
            mpmath variable.

        """

        if not(isinstance(initial_value, numpy.ndarray)):
            initial_value = numpy.array(initial_value)

        if initial_value.size == 1:
            mpc = mpmath.mpc(complex(initial_value.real) + 1j * complex(initial_value.imag))
        else:
            mpc = mpmath.mpc(complex(1) + 1j * complex(0)) * initial_value

        return mpc

    def exponentiate(self, power):
        """Exponentiates to the power.

        Parameters
        ----------
        power :
            The power to raise to.

        Returns
        -------
        type
            Exponential.

        """
        return mpmath_exp(power)

    def sin(self, power):
        """Sin to the power.

        Parameters
        ----------
        power :
            The power to raise to.

        Returns
        -------
        type
            Sin.

        """
        return mpmath_sin(power)

    def cos(self, power):
        """Cos to the power.

        Parameters
        ----------
        power :
            The power to raise to.

        Returns
        -------
        type
            Cos.

        """
        return mpmath_cos(power)

    def toComplex(self, variable):
        """Converts calculation variable to native python complex.

        Parameters
        ----------
        variable :
            Calculation variable to convert.

        Returns
        -------
        type
            Native python complex variable.

        """
        # return complex(variable)

        return numpy.array(variable, dtype=complex)


class CalculationStrategyMath(CalculationStrategy):
    """Use plain python for calculation."""
    def createVariable(self, initial_value):
        """Factory method for calculation variable.

        Parameters
        ----------
        initial_value :
            Initial value of the variable.

        Returns
        -------
        type
            mpmath variable.

        """
        return initial_value + 0j # complex(initial_value)

    def exponentiate(self, power):
        """Exponentiates to the power.

        Parameters
        ----------
        power :
            The power to raise to.

        Returns
        -------
        type
            Exponential.

        """
        try:
            ans =  numpy.exp(power)
        except:
            ans = float("Inf")
        return ans

    def sin(self, power):
        """Sin to the power.

        Parameters
        ----------
        power :
            The power to raise to.

        Returns
        -------
        type
            Sin.

        """
        return numpy.sin(power)

    def cos(self, power):
        """Cos to the power.

        Parameters
        ----------
        power :
            The power to raise to.

        Returns
        -------
        type
            Cos.

        """
        return numpy.cos(power)


    def toComplex(self, variable):
        """Converts calculation variable to native python complex.

        Parameters
        ----------
        variable :
            Calculation variable to convert.

        Returns
        -------
        type
            Native python complex variable.

        """
        return complex(variable)


class PerfectCrystalDiffraction(object):
    """
    PerfectCrystalDiffraction is the calculator of the perfect crystal.

    Two steps:
    * create the PerfectCrystalDiffraction instance with the crystal data (usually picked up from DiffractionSetup).
    * call PerfectCrystalDiffraction.calculateDiffraction( photon(s) )

    Notes:
    * arrays can be used, but compatible arrays in photons and in bragg_angle, psi_0, psi_H, psi_H, etc.
    * Both sigma and pi amplitudes can be calculated in the same call, not need to create different
      instances for sigma and pi.

    Constructor.

    Parameters
    ----------
    geometry_type: instance of BraggDiffraction, LaueDiffraction, BraggTransmission, or LaueTransmission
    bragg_normal : instance of Vector
        The H vector.
    surface_normal : instance of Vector
        The n vector.
    bragg_angle : float or numpy array
        The Bragg angle(s).
    psi_0 : complex of numpy array
        The structire factor Psi0.
    psi_H: complex of numpy array
        The structire factor PsiH.
    psi_H_bar: complex of numpy array
        The structire factor Psi(-H).
    thickness : float
        the crystal thickness in m.
    d_spacing : float
        the crystal dSpacing in m.

    """
    isDebug = False

    def __init__(self, geometry_type, bragg_normal, surface_normal, bragg_angle, psi_0, psi_H, psi_H_bar, thickness, d_spacing):

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

    #
    # getters
    #
    def braggNormal(self):
        """Returns the Bragg normal, i.e. normal on the reflection planes with modulus 2 pi / d_spacing.
        :return: Bragg normal.

        Parameters
        ----------

        Returns
        -------

        """
        return self._bragg_normal

    def surfaceNormal(self):
        """Returns the surface normal that points outwards the crystal.
        :return: Surface normal.

        Parameters
        ----------

        Returns
        -------

        """
        return self._surface_normal

    def surfaceNormalInwards(self):
        """Returns the surface normal that points inwards the crystal.
        :return: Surface normal.

        Parameters
        ----------

        Returns
        -------

        """
        return self._surface_normal.scalarMultiplication(-1.0)

    def braggAngle(self):
        """Returns the Bragg angle.
        :return: The Bragg angle.

        Parameters
        ----------

        Returns
        -------

        """
        return self._bragg_angle

    def Psi0(self):
        """Returns Psi0 as defined in Zachariasen [3-95].
        :return: Psi0.

        Parameters
        ----------

        Returns
        -------

        """
        return self._psi_0

    def PsiH(self):
        """Returns Psi0 as defined in Zachariasen [3-95].
        :return: PsiH.

        Parameters
        ----------

        Returns
        -------

        """
        return self._psi_H

    def PsiHBar(self):
        """Returns Psi0 as defined in Zachariasen [3-95].
        :return: PsiHBar.

        Parameters
        ----------

        Returns
        -------

        """
        return self._psi_H_bar

    def thickness(self):
        """Returns crystal thickness.
        :return: Thickness of the crystal.

        Parameters
        ----------

        Returns
        -------

        """
        return self._thickness

    def dSpacing(self):
        """Returns distance between the reflection planes.
        :return: Distance between the reflection planes.

        Parameters
        ----------

        Returns
        -------

        """
        return self._d_spacing

    def geometryType(self):
        """Returns the geometry types, i.e. BraggTransmission, LaueDiffraction,...
        :return: Geometry type.

        Parameters
        ----------

        Returns
        -------

        """
        return self._geometryType

    #
    # i/o
    #
    def log(self, string):
        """Logs a string.

        Parameters
        ----------
        string :
            String to log.

        Returns
        -------

        """
        print(string)

    def logDebug(self, string):
        """Logs a debug string.

        Parameters
        ----------
        string :
            String to log.

        Returns
        -------

        """
        self.log("<DEBUG>: " + string)

    def _logMembers(self, zac_b, zac_alpha, photon_in, photon_out, result):
        """Debug logs the member variables and other relevant partial results.

        Parameters
        ----------
        zac_b :
            Asymmetry ratio b
        zac_alpha :
            Diffraction index difference of crystal fields.
        photon_in :
            Incoming photon.
        result :
            Resulting complex amplitudes of the diffraction/transmission.
        photon_out :


        Returns
        -------

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
        self.logDebug('intensity S: ' + str(numpy.abs(result["S"]) ** 2))
        self.logDebug('intensity P: ' + str(numpy.abs(result["P"]) ** 2))

    #
    # basic auxiliar calculations
    #
    def _calculateGamma(self, photon):
        """Calculates the projection cosine gamma as defined in Zachariasen [3-115].

        Parameters
        ----------
        photon :
            Photon that is projected onto the surface normal.

        Returns
        -------
        type
            Projection cosine gamma.

        """
        # gamma = photon.unitDirectionVector().scalarProduct(self.surfaceNormal().getNormalizedVector())
        # # Our crystal normal is pointing outside the crystal medium. Zachariasen's normal points
        # # into the crystal medium (pag 112). Therefore, we change the sign.
        # gamma = -gamma
        # return gamma

        return photon.unitDirectionVector().scalarProduct(self.surfaceNormalInwards())

    def _calculatePhotonOut(self, photon_in):
        """Solves the scattering equation to calculates the outgoing photon from the incoming photon and the Bragg normal.

        1) Calculates the parammemm component of K: k_out_par = k_in_par + H_par
        2) Uses the conservation of the wavevector modulus to calculate the k_out_perp

        It is valid for diffraction not at the Bragg angle.

        Parameters
        ----------
        photon_in :
            Incoming photon or photon bunch.

        Returns
        -------
        instance of ComplexAmplitudePhoton
            Outgoing photon or photon bunch

        """
        # # Retrieve k_0.
        # k_in = photon_in.wavevector()

        # # Solve unscaled Laue equation.
        # k_out = self.braggNormal().addVector(k_in)

        # Create photon in k_out direction and scale by setting the photon energy.
        # photon_out = Photon(photon_in.energy(), k_out)
        # """
        # GENERAL VERSION:
        # Solves the Laue equation for the parallel components of the vectors and
        # uses the conservation of the wavevector modulus to calculate the outgoing wavevector
        # even for diffraction not at the Bragg angle.
        # """
        # Retrieve k_0.
        k_in = photon_in.wavevector()

        # Decompose the vector into a component parallel to the surface normal and
        # a component parallel to the surface: (k_in * n) n.
        k_in_normal = self.surfaceNormal().scalarMultiplication(k_in.scalarProduct(self.surfaceNormal()))
        k_in_parallel = k_in.subtractVector(k_in_normal)

        # Retrieve the B_H vector.
        B_H = self.braggNormal()

        # Decompose the vector into a component parallel to the surface normal and
        # a component parallel to the surface: (B_H * n) n.
        B_H_normal = self.surfaceNormal().scalarMultiplication(B_H.scalarProduct(self.surfaceNormal()))
        B_H_parallel = B_H.subtractVector(B_H_normal)

        # Apply the Laue formula for the parallel components.
        k_out_parallel = k_in_parallel.addVector(B_H_parallel)

        # Calculate K_out normal.
        k_out_normal_modulus = numpy.sqrt(k_in.norm() ** 2 - k_out_parallel.norm() ** 2)
        k_out_normal = self.surfaceNormal().scalarMultiplication(k_out_normal_modulus)

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
            self.logDebug("surface normal" + str(self.surfaceNormal().components()))
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
        """Calculates alpha ("refraction index difference between waves in the crystal") as defined in Zachariasen [3-114b].

        Parameters
        ----------
        photon_in :
            Incoming photon.

        Returns
        -------
        type
            alpha.

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
        """Calculates alpha ("refraction index difference between waves in the crystal") as defined in Zachariasen [3-114b].

        Parameters
        ----------
        photon_in :
            Incoming photon.

        Returns
        -------
        type
            alpha.

        """
        k0_dot_H = photon_in.wavevector().scalarProduct(self.braggNormal()) # scalar product k0 and H.
        wavenumber = photon_in.wavenumber() #  norm of k0.
        alpha = - (wavenumber ** -2) * (self.braggNormal().norm() ** 2 + 2 * k0_dot_H)
        return alpha


    def _calculateZacB(self, photon_in, photon_out):
        """Calculates asymmetry ratio b as defined in Zachariasen [3-115].

        Parameters
        ----------
        photon_in :
            Incoming photon.
        photon_out :
            Outgoing photon.

        Returns
        -------
        type
            Asymmetry ratio b.

        """
        # TODO: revise this algorithm, it is not exactly as in Zachariasen [3-115]
        numerator   = self.surfaceNormal().scalarProduct(photon_in.wavevector())
        denominator = self.surfaceNormal().scalarProduct(photon_out.wavevector())
        zac_b = numerator / denominator

        return zac_b


    def _calculateGuigayB(self, photon_in):
        """Calculates asymmetry ratio b as defined in Guigay.

        Parameters
        ----------
        photon_in :
            Incoming photon.

        Returns
        -------
        type
            Asymmetry ratio b.
            Note that this b changes when K0 (photon_in.wavevector()) changes

        """
        KH = photon_in.wavevector().addVector(self.braggNormal())
        photon_outG = Photon(energy_in_ev=photon_in.energy(), direction_vector=KH)
        return self._calculateGamma(photon_in) / self._calculateGamma(photon_outG)


    def _calculateZacQ(self, zac_b, effective_psi_h, effective_psi_h_bar):
        """Calculates q as defined in Zachariasen [3-123].

        Parameters
        ----------
        zac_b :
            Asymmetry ratio b as defined in Zachariasen [3-115].
        effective_psi_h :
            Effective PsiH (depending of polarisation. See text following [3.-139]).
        effective_psi_h_bar :
            Effective PsiHBar (depending of polarisation. See text following [3.-139]).

        Returns
        -------
        type
            q.

        """
        return zac_b * effective_psi_h * effective_psi_h_bar

    def _calculateZacZ(self, zac_b, zac_alpha):
        """Calcualtes z as defined in Zachariasen [3-123].

        Parameters
        ----------
        zac_b :
            Asymmetry ratio b as defined in Zachariasen [3-115].
        zac_alpha :
            Diffraction index difference of crystal fields.

        Returns
        -------
        type
            z.

        """
        return (1.0e0 - zac_b) * 0.5e0 * self.Psi0() + zac_b * 0.5e0 * zac_alpha

    #
    # math tools
    #
    def _createVariable(self, initial_value):
        """Factory method for calculation variable. Delegates to active calculation strategy.

        Parameters
        ----------
        initial_value :
            Inital value of the variable.

        Returns
        -------
        type
            Variable to use for the calculation.

        """
        return self._calculation_strategy.createVariable(initial_value)

    def _exponentiate(self, power):
        """Exponentiates to the power using active calculation strategy. (plain python or arbitrary precision)

        Parameters
        ----------
        power :
            Calculation variable.

        Returns
        -------
        type
            Exponential.

        """
        return self._calculation_strategy.exponentiate(self._createVariable(power))

    def _sin(self, power):
        """Sin to the power using active calculation strategy. (plain python or arbitrary precision)

        Parameters
        ----------
        power :
            Calculation variable.

        Returns
        -------
        type
            Sin.

        """
        return self._calculation_strategy.sin(self._createVariable(power))

    def _cos(self, power):
        """Cos to the power using active calculation strategy. (plain python or arbitrary precision)

        Parameters
        ----------
        power :
            Calculation variable.

        Returns
        -------
        type
            Cos.

        """
        return self._calculation_strategy.cos(self._createVariable(power))

    def _toComplex(self, variable):
        """Converts calculation variable to complex. Delegates to active calculation strategy.

        Parameters
        ----------
        variable :
            Calculation variable.

        Returns
        -------
        type
            Calculation variable as complex.

        """
        return self._calculation_strategy.toComplex(variable)


    #
    # final auxiliar parameters calculations
    #
    def _calculateComplexAmplitude(self, photon_in, zac_q, zac_z, gamma_0, effective_psi_h_bar):
        """Calculates the complex amplitude of the questioned wave: diffracted or transmission.

        Parameters
        ----------
        photon_in :
            Incoming photon.
        zac_q :
            q as defined in Zachariasen [3-123].
        zac_z :
            z as defined in Zachariasen [3-123].
        gamma_0 :
            Projection cosine as defined in Zachariasen [3-115].
        effective_psi_h_bar :
            Effective PsiHBar (depending of polarisation. See text following [3.-139]).

        Returns
        -------
        type
            Complex amplitude.

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
        """Calculates complex amplitude for the S polarization.

        Parameters
        ----------
        photon_in :
            Incoming photon.
        zac_z :
            z as defined in Zachariasen [3-123].
        gamma_0 :
            Projection cosine as defined in Zachariasen [3-115].
        zac_b :
            

        Returns
        -------
        type
            Complex amplitude of S polarization.

        """
        zac_q = self._calculateZacQ(zac_b, self.PsiH(), self.PsiHBar())
        return self._calculateComplexAmplitude(photon_in, zac_q, zac_z, gamma_0, self.PsiHBar())

    def _calculatePolarizationP(self, photon_in, zac_b, zac_z, gamma_0):
        """Calculates complex amplitude for the P polarization.

        Parameters
        ----------
        photon_in :
            Incoming photon.
        zac_b :
            Asymmetry ratio b as defined in Zachariasen [3-115].
        zac_z :
            z as defined in Zachariasen [3-123].
        gamma_0 :
            Projection cosine as defined in Zachariasen [3-115].

        Returns
        -------
        type
            Complex amplitude of P polarization.

        """
        effective_psi_h = self.PsiH() * numpy.cos(2 * self.braggAngle())
        effective_psi_h_bar = self.PsiHBar() * numpy.cos(2 * self.braggAngle())

        zac_q = self._calculateZacQ(zac_b, effective_psi_h, effective_psi_h_bar)
        return self._calculateComplexAmplitude(photon_in, zac_q, zac_z, gamma_0, effective_psi_h_bar)

    #
    # final parameters calculations
    #
    def calculateDiffraction(self,
                             photon_in,
                             calculation_method=0, # 0=Zachariasen, 1=Guigay
                             is_thick=0, # for Guigay only
                             use_transfer_matrix=0, # for Guigay only
                             ):
        """Calculate diffraction for incoming photon.

        Parameters
        ----------
        photon_in : instance of ComplexAmplitudePhoton
            Incoming photon or Photon bunch.
        calculation_method : int
            0 : Zachariasen, 1 : Guigay
        is_thick : int
            0=No, 1=Yes (for calculation_method=1 only)
        use_transfer_matrix : int
            0=No, 1=Yes (for calculation_method=1 only)

        Returns
        -------
        dict
            The complex amplitudes of the diffraction weighted for power for the two polarizations are found
            in the kwys ["S"]  and ["P"].

        """

        if calculation_method == 0:
            # print(">>>> Using Zachariasen equations...")
            return self.calculateDiffractionZachariasen(photon_in)
        else:
            # print(">>>> Using Guigay equations...")
            return self.calculateDiffractionGuigay(photon_in, is_thick=is_thick, use_transfer_matrix=use_transfer_matrix)


    def calculateDiffractionZachariasen(self, photon_in):
        """Calculate diffraction for incoming photon.

        Parameters
        ----------
        photon_in :
            Incoming photon.

        Returns
        -------
        dict
            The complex amplitudes of the diffraction weighted for power for the two polarizations are found
            in the kwys ["S"]  and ["P"].

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
        """Calculate diffraction for incoming photon.

        Parameters
        ----------
        photon_in : instance of ComplexAmplitudePhoton
            Incoming photon.
        debug : int
             0=No, 1=Yes.
        s_ratio : float
             the sin(theta)/lambda ratio (Default value = calculated from Bragg law)
        is_thick : int
             0=No, 1=Yes.
        use_transfer_matrix : int
             0=No, 1=Yes (Default value = 0). It is is faster to use use_transfer_matrix=0

        Returns
        -------
        dict
            The complex amplitudes of the diffraction weighted for power for the two polarizations are found
            in the kwys ["S"]  and ["P"].
            Other optional parameters are found at the keys:
            * sigma-polarized reflectivity: "s"
            * pi-polarized reflectivity: "p"
            * Transfer matrix: "m11_s" "m12_s" "m21_s" "m22_s" "m11_p" "m12_p" "m21_p" "m22_p"
            * Scattering matrix: "s11_s" "s12_s" "s21_s" "s22_s" "s11_p" "s12_p" "s21_p" "s22_p"
            * "gamma_0"
            * "alpha"
            * "b"

        """
        # Initialize return variable.

        result = {"S": None,
                  "P": None}

        guigay_b = self._calculateGuigayB(photon_in)  # gamma_0 / gamma_H
        if debug: print("guigay_b: ", guigay_b)

        if use_transfer_matrix:
            transfer_matrix_s = self.calculateTransferMatrix(photon_in, polarization=0, is_thick=is_thick)
            m11_s, m12_s, m21_s, m22_s = transfer_matrix_s
            scattering_matrix_s = self.calculateScatteringMatrixFromTransferMatrix(transfer_matrix_s)
            s11_s, s12_s, s21_s, s22_s = scattering_matrix_s

            transfer_matrix_p = self.calculateTransferMatrix(photon_in, polarization=1, is_thick=is_thick)
            m11_p, m12_p, m21_p, m22_p = transfer_matrix_p
            scattering_matrix_p = self.calculateScatteringMatrixFromTransferMatrix(transfer_matrix_p)
            s11_p, s12_p, s21_p, s22_p = scattering_matrix_p

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

        else:
            # Calculate crystal field refraction index difference.
            # Note that Guigay's definition of alpha has the opposite sign as in Zachariasen!
            alpha = self._calculateGuigayAlpha(photon_in)
            if debug: print("guigay alpha: ", alpha)

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



        # Calculate complex amplitude for S and P polarization.


        # store results

        if (self.geometryType() == BraggDiffraction() or self.geometryType() == LaueDiffraction()):
            result["S"] = complex_amplitude_s / numpy.sqrt(abs(guigay_b))
            result["P"] = complex_amplitude_p / numpy.sqrt(abs(guigay_b))
        else:
            result["S"] = complex_amplitude_s
            result["P"] = complex_amplitude_p

        if use_transfer_matrix:
            result["s"] = complex_amplitude_s
            result["p"] = complex_amplitude_p
            result["alpha"] = alpha
            result["b"] = guigay_b

        # If debugging output is turned on.
        if self.isDebug:
            self._logMembers(guigay_b, alpha, photon_in, self._calculatePhotonOut(photon_in), result)

        # Returns the complex amplitudes.
        return result

    def calculateTransferMatrix(self, photon_in, polarization=0, is_thick=0):
        """

        Parameters
        ----------
        photon_in :
            
        polarization :
             (Default value = 0)
        is_thick :
             (Default value = 0)

        Returns
        -------

        """

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
        """

        Parameters
        ----------
        transfer_matrix :
            

        Returns
        -------

        """
        m11, m12, m21, m22 = transfer_matrix
        s11 = m11 - m12 * m21 / m22
        s12 = m12 / m22
        s21 = -m21 / m22
        s22 = 1 / m22
        return s11, s12, s21, s22

    def calculateScatteringMatrix(self, photon_in, polarization=0):
        """

        Parameters
        ----------
        photon_in :
            
        polarization :
             (Default value = 0)

        Returns
        -------

        """
        transfer_matrix = self.calculateTransferMatrix(photon_in, polarization=polarization)
        return self.calculateScatteringMatrixFromTransferMatrix(transfer_matrix)


if __name__ == "__main__":
    a = CalculationStrategyMPMath()

    a0 = a.createVariable(numpy.array(0))
    print("cos(0)", a.cos(a0))

    api = a.createVariable(numpy.array(numpy.pi))
    print("cos(pi)", a.cos(api))

    api = a.createVariable(numpy.array([numpy.pi] * 10))
    print("cos(pi)", a.cos(api))
    # print("exp(pi)", a.exponentiate(api))

