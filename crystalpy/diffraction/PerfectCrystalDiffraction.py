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

        Raises
        ------
        Exception
            Must override this method.

        """
        raise Exception("Must override this method.")

    def exponentiate(self, power):
        """Exponentiates to the power.

        Parameters
        ----------
        power :
            The power to raise to.

        Raises
        ------
        Exception
            Must override this method.

        """
        raise Exception("Must override this method.")

    def sin(self, power):
        """Sin to the power.

        Parameters
        ----------
        power :
            The power to raise to.

        Raises
        ------
        Exception
            Must override this method.

        """
        raise Exception("Must override this method.")

    def cos(self, power):
        """Cos to the power.

        Parameters
        ----------
        power :
            The power to raise to.

        Raises
        ------
        Exception
            Must override this method.

        """
        raise Exception("Must override this method.")


    def toComplex(self, variable):
        """Converts calculation variable to native python complex.

        Parameters
        ----------
        variable :
            Calculation variable to convert.

        Raises
        ------
        Exception
            Must override this method.

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
        initial_value : float, complex or numpy array
            Initial value of the variable.

        Returns
        -------
        instance of CalculationStrategyMPMath
            variable.

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
        power : float
            The power to raise to.

        Returns
        -------
        mpmath variable
            Exponential.

        """
        return mpmath_exp(power)

    def sin(self, power):
        """Sin function.

        Parameters
        ----------
        power : float or numpy array
            The arg of sin.

        Returns
        -------
        mpmath variable
            Sin.

        """
        return mpmath_sin(power)

    def cos(self, power):
        """Cos function.

        Parameters
        ----------
        power : float or numpy array
            The arg of cos.

        Returns
        -------
        mpmath variable
            Cos.

        """
        return mpmath_cos(power)

    def toComplex(self, variable):
        """Converts calculation variable to native python complex.

        Parameters
        ----------
        variable :
            variable to convert.

        Returns
        -------
        numpy array
            Native python complex variable.

        """
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
        instance of CalculationStrategyMath
            variable.

        """
        return initial_value + 0j # complex(initial_value)

    def exponentiate(self, power):
        """Exponentiates to the power.

        Parameters
        ----------
        power : float
            The power to raise to.

        Returns
        -------
        numpy array
            Exponential.

        """
        try:
            ans =  numpy.exp(power)
        except:
            ans = float("Inf")
        return ans

    def sin(self, power):
        """Sin function.

        Parameters
        ----------
        power :
            The sin argument.

        Returns
        -------
        numpy array
            Sin.

        """
        return numpy.sin(power)

    def cos(self, power):
        """Cos function.

        Parameters
        ----------
        power :
            The coa argument.

        Returns
        -------
        numpy array
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
        numpy array (complex)
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


    @classmethod
    def initializeFromDiffractionSetupAndEnergy(cls, diffraction_setup, energy,
                                 geometry_type=None,
                                 bragg_normal=None,
                                 surface_normal=None,
                                 # bragg_angle=None,
                                 # psi_0=None,
                                 # psi_H=None,
                                 # psi_H_bar=None,
                                 thickness=None,
                                 d_spacing=None,
                                 ):
        """
        Creates a PerfectCrystalDiffraction instance from parameters in a DiffractionSetupAbstract instance and a
        photon energy array.

        Parameters
        ----------
        diffraction_setup : instance of PerfectCrystalDiffraction

        energy : numpy array

        geometry_type: instance of BraggDiffraction, LaueDiffraction, BraggTransmission, or LaueTransmission

        bragg_normal : instance of Vector, optional
            if None, retrieve from DiffractionSetup

        surface_normal : instance of Vector, optional
            if None, retrieve from DiffractionSetup

        thickness : float or numpy array, optional
            crystal thickness in m. If None, retrieve from DiffractionSetup

        d_spacing : float or numpy array
            d-spacing in m. If None, retrieve from DiffractionSetup


        Returns
        -------
        PerfectCrystalDiffraction instance


        """

        #
        # energy-depending variables
        #

        # Retrieve bragg angle.
        bragg_angle = diffraction_setup.angleBragg(energy)

        psi_0, psi_H, psi_H_bar = diffraction_setup.psiAll(energy)

        geometry_type = geometry_type if not geometry_type is None else diffraction_setup.geometryType()
        bragg_normal = bragg_normal if not bragg_normal is None else diffraction_setup.vectorH()
        surface_normal = surface_normal if not surface_normal is None else diffraction_setup.vectorNormalSurface()
        thickness = thickness if not thickness is None else diffraction_setup.thickness()
        d_spacing = d_spacing if not d_spacing is None else diffraction_setup.dSpacingSI()

        return PerfectCrystalDiffraction(
            geometry_type,
            bragg_normal,
            surface_normal,
            bragg_angle,
            psi_0,
            psi_H,
            psi_H_bar,
            thickness,
            d_spacing)


    #
    # getters
    #
    def braggNormal(self):
        """Returns the Bragg normal, i.e. normal on the reflection planes with modulus 2 pi / d_spacing.

        Returns
        -------
        instance of Vector
            The H vector (normal to Bragg planes, modulus 2 pi / d_spacing in m^-1))

        """
        return self._bragg_normal

    def surfaceNormal(self):
        """Returns the surface normal that points outwards the crystal.

        Returns
        -------
        instance of Vector
            The normal to the surface.

        """
        return self._surface_normal

    def surfaceNormalInwards(self):
        """Returns the surface normal that points inwards the crystal.

        Returns
        -------
        instance of Vector
            The normal to the surface.

        """
        return self._surface_normal.scalarMultiplication(-1.0)

    def braggAngle(self):
        """Returns the Bragg angle.

        Returns
        -------
        float or numpy array
            The Bragg angle in rad.

        """
        return self._bragg_angle

    def Psi0(self):
        """Returns Psi0 as defined in Zachariasen [3-95].

        Returns
        -------
        complex or numpy array
            The Psi0 value.

        """
        return self._psi_0

    def PsiH(self):
        """Returns PsiH as defined in Zachariasen [3-95].

        Parameters
        ----------

        Returns
        -------
        complex or numpy array
            The PsiH value.

        """
        return self._psi_H

    def PsiHBar(self):
        """Returns PsiHBar as defined in Zachariasen [3-95].

        Returns
        -------
        complex or numpy array
            The PsiHBar value.

        """
        return self._psi_H_bar

    def thickness(self):
        """Returns crystal thickness.

        Returns
        -------
        float or numpy array
            The thickness of the crystal in m.

        """
        return self._thickness

    def dSpacing(self):
        """Returns the distance between the reflection planes in A.

        Returns
        -------
        float or numpy array
            Distance between the reflection planes in A.

        """
        return self._d_spacing

    def geometryType(self):
        """Returns the geometry types, i.e. BraggTransmission, LaueDiffraction,...

        Returns
        -------
                geometry_type: instance of BraggDiffraction, LaueDiffraction, BraggTransmission, or LaueTransmission

        """
        return self._geometryType

    #
    # i/o
    #
    def log(self, string):
        """Logs (prints) a string.

        Parameters
        ----------
        string :
            String to log.

        """
        print(string)

    def logDebug(self, string):
        """Logs (prints) a debug string.

        Parameters
        ----------
        string :
            String to log.

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
        photon : instance of ComplexAmplitudePhoton
            Photon or photon bunch that is projected onto the surface normal.

        Returns
        -------
        float or numpy array
            Projection cosine gamma.

        """
        return photon.unitDirectionVector().scalarProduct(self.surfaceNormalInwards())

    def _calculatePhotonOut(self, photon_in,
                            method=1,
                            apply_reflectivity=False,
                            calculation_method=0,  # 0=Zachariasen, 1=Guigay
                            is_thick=0,  # for Guigay only
                            use_transfer_matrix=0,  # for Guigay only
                            ):
        """
        Solves the scattering equation to calculates the outgoing photon from the incoming photon and the Bragg normal.

        In case of diffracted photon (Laue or Bragg)
        1) Calculates the parallel component of K: k_out_par = k_in_par + H_par
        2) Uses the conservation of the wavevector modulus to calculate the k_out_perp

        todo: In case of forward diffracted (transmitted) Bragg or Laue:

        It is valid for diffraction not at the Bragg angle.

        Warnings
        --------
        Note that the new photon has correct direction, but the complex amplitudes are not changed, just
        copied from the photon in.

        Parameters
        ----------
        photon_in : instance of ComplexAmplitudePhoton
            Incoming photon or photon bunch.
        method : int
            select the calculated method: 0=old (to be deleted), 1=new

        Returns
        -------
        instance of ComplexAmplitudePhoton
            Outgoing photon or photon bunch

        """
        # GENERAL VERSION:
        # Solves the Laue equation for the parallel components of the vectors and
        # uses the conservation of the wavevector modulus to calculate the outgoing wavevector
        # even for diffraction not at the Bragg angle.

        # Retrieve k_0.
        k_in = photon_in.wavevector()
        # Retrieve the B_H vector.
        B_H = self.braggNormal()

        if method == 0: # old method, todo: delete

            # Decompose the vector into a component parallel to the surface normal and
            # a component parallel to the surface: (k_in * n) n.
            k_in_normal = self.surfaceNormal().scalarMultiplication(k_in.scalarProduct(self.surfaceNormal()))
            k_in_parallel = k_in.subtractVector(k_in_normal)


            # Decompose the vector into a component parallel to the surface normal and
            # a component parallel to the surface: (B_H * n) n.
            B_H_normal = self.surfaceNormal().scalarMultiplication(B_H.scalarProduct(self.surfaceNormal()))
            B_H_parallel = B_H.subtractVector(B_H_normal)

            # Apply the Laue formula for the parallel components.
            k_out_parallel = k_in_parallel.addVector(B_H_parallel)

            # Calculate K_out normal.
            k_out_normal_modulus = numpy.sqrt(k_in.norm() ** 2 - k_out_parallel.norm() ** 2)
            k_out_normal = self.surfaceNormal().scalarMultiplication(k_out_normal_modulus)

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

        elif method == 1: # new method
            if self.geometryType() == BraggDiffraction():
                k_out = k_in.scatteringOnSurface(self.surfaceNormal(), B_H, use_sign_of=+1)
            elif self.geometryType() == LaueDiffraction():
                k_out = k_in.scatteringOnSurface(self.surfaceNormal(), B_H, use_sign_of=-1)
            elif self.geometryType() == BraggTransmission():
                k_out = k_in.scatteringOnSurface(self.surfaceNormal(), B_H, use_sign_of=+1) # todo: fix
            elif self.geometryType() == LaueTransmission():
                k_out = k_in.scatteringOnSurface(self.surfaceNormal(), B_H, use_sign_of=-1) # todo: fix


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


        if apply_reflectivity:
            coeffs = self.calculateDiffraction(photon_in,
                                                          calculation_method=calculation_method,
                                                          is_thick=is_thick,
                                                          use_transfer_matrix=use_transfer_matrix)

            # apply reflectivities
            photon_out.rescaleEsigma(coeffs["S"])
            photon_out.rescaleEpi(coeffs["P"])

        # Return outgoing photon.

        return photon_out


    def _calculateAlphaZac(self, photon_in):
        """
        Calculates alpha ("refraction index difference between waves in the crystal") as defined in Zachariasen [3-114b].

        Parameters
        ----------
        photon_in : instance of ComplexAmplitudePhoton
            Incoming photon or photon bunch.

        Returns
        -------
        float or numpy array
            The alpha value.

        """
        # Calculate scalar product k_0 and B_H.
        k_0_times_B_h = photon_in.wavevector().scalarProduct(self.braggNormal())

        # Get norm k_0.
        wavenumber = photon_in.wavenumber()

        # Calculate alpha.
        zac_alpha = (wavenumber ** -2) * (self.braggNormal().norm() ** 2 + 2 * k_0_times_B_h)

        # Return alpha.
        return zac_alpha

    def _calculateAlphaGuigay(self, photon_in):
        """Calculates alpha ("refraction index difference between waves in the crystal") as defined in Guigay eq. XXX.

        It is the same as the Zachariasen alpha value with the opposite sign.

        Parameters
        ----------
        photon_in : instance of ComplexAmplitudePhoton
            Incoming photon or photon bunch.

        Returns
        -------
        float or numpy array
            The alpha value.

        """
        k0_dot_H = photon_in.wavevector().scalarProduct(self.braggNormal()) # scalar product k0 and H.
        wavenumber = photon_in.wavenumber() #  norm of k0.
        alpha = - (wavenumber ** -2) * (self.braggNormal().norm() ** 2 + 2 * k0_dot_H)
        return alpha


    def _calculateZacB(self, photon_in, photon_out):
        """Calculates asymmetry ratio b as defined in Zachariasen equation [3-114a].

        Parameters
        ----------
        photon_in : instance of ComplexAmplitudePhoton
            Incoming photon or photon bunch.
        photon_out :
            Outgoing photon.

        Returns
        -------
        float or numpy array
            Asymmetry ratio b.

        """
        numerator   = self.surfaceNormalInwards().scalarProduct(self.braggNormal())
        denominator = self.surfaceNormalInwards().scalarProduct(photon_in.wavevector())
        zac_b = 1.0 / (numerator / denominator + 1)
        return zac_b


    def _calculateGuigayB(self, photon_in):
        """Calculates asymmetry ratio b as defined in Guigay equation 17.

        Parameters
        ----------
        photon_in : instance of ComplexAmplitudePhoton
            Incoming photon or photon bunch.

        Returns
        -------
        float or numpy array
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
        complex or numpy array
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
        complex or numpy array
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
        photon_in : instance of ComplexAmplitudePhoton
            Incoming photon or photon bunch.
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
        photon_in : instance of ComplexAmplitudePhoton
            Incoming photon or photon bunch.
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
        photon_in : instance of ComplexAmplitudePhoton
            Incoming photon or photon bunch.
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
        photon_in : instance of ComplexAmplitudePhoton
            Incoming photon or photon bunch.

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
        zac_b = self._calculateZacB(photon_in, photon_out)
        # zac_b = self._calculateGuigayB(photon_in)  # todo: check if this is the same for Zac

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
            result["S"] *= (1.0 / numpy.sqrt(abs(zac_b)))
            result["P"] *= (1.0 / numpy.sqrt(abs(zac_b)))

        # If debugging output is turned on.
        if self.isDebug:
            self._logMembers(zac_b, zac_alpha, photon_in, photon_out, result)

        # Returns the complex amplitudes.
        return result

    def calculateDiffractionGuigay(self, photon_in,
                                   debug=0,
                                   s_ratio=None,
                                   is_thick=0,
                                   use_transfer_matrix=0, # is faster to use use_transfer_matrix=0
                                   ):
        """Calculate diffraction for incoming photon.

        Parameters
        ----------
        photon_in : instance of ComplexAmplitudePhoton
            Incoming photon or photon bunch.
        debug : int
             0=No, 1=Yes.
        s_ratio : float
             the sin(theta)/lambda ratio (if None (default) it is calculated from Bragg law).
        is_thick : int
             0=No, 1=Yes.
        use_transfer_matrix : int
             0=No, 1=Yes. It is is faster to use use_transfer_matrix=0.

        Returns
        -------
        dict
            The complex amplitudes of the diffraction weighted for power for the two polarizations are found
            in the kwys ["S"]  and ["P"].
            If use_transfer_matrix=1, other optional parameters are found at the keys:
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
        alpha = self._calculateAlphaGuigay(photon_in)
        gamma_0 = self._calculateGamma(photon_in)
        T = self.thickness() / gamma_0
        if debug:
            print("guigay_b: ", guigay_b)
            print("guigay alpha: ", alpha)
            print("gamma_0: ", gamma_0)
            print("T: ", T)

        if use_transfer_matrix:
            transfer_matrix_s = self.calculateTransferMatrix(photon_in, polarization=0, is_thick=is_thick,
                                                             alpha=alpha, guigay_b=guigay_b, T=T)
            m11_s, m12_s, m21_s, m22_s = transfer_matrix_s
            scattering_matrix_s = self.calculateScatteringMatrixFromTransferMatrix(transfer_matrix_s)
            s11_s, s12_s, s21_s, s22_s = scattering_matrix_s

            transfer_matrix_p = self.calculateTransferMatrix(photon_in, polarization=1, is_thick=is_thick,
                                                             alpha=alpha, guigay_b=guigay_b, T=T)
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

            result["s"] = complex_amplitude_s
            result["p"] = complex_amplitude_p
            result["alpha"] = alpha
            result["b"] = guigay_b

        else:

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



        # store complex amplitude weighted for power for S and P polarization.

        if (self.geometryType() == BraggDiffraction() or self.geometryType() == LaueDiffraction()):
            result["S"] = complex_amplitude_s / numpy.sqrt(abs(guigay_b))
            result["P"] = complex_amplitude_p / numpy.sqrt(abs(guigay_b))
        else:
            result["S"] = complex_amplitude_s
            result["P"] = complex_amplitude_p

        # If debugging output is turned on.
        if self.isDebug:
            self._logMembers(guigay_b, alpha, photon_in, self._calculatePhotonOut(photon_in), result)

        # Returns the complex amplitudes.
        return result

    def calculateTransferMatrix(self, photon_in, polarization=0, is_thick=0, alpha=None, guigay_b=None, T=None):
        """
        Calculates the transfer matrix (see equation XXX in Guigay and Sanchez del Rio).

        Parameters
        ----------
        photon_in : instance of ComplexAmplitudePhoton
            Incoming photon or photon bunch.
        polarization : int
             0=sigma, 1=pi.
        is_thick : int
             Use thick crystal approximation: 0=No, 1=Yes.
        alpha : float or numpy array
            The alpha value (if None, it is internally calculated).
        guigay_b : float or numpy array
            The b (asymmetry factor) value (if None, it is internally calculated).
        T : float or numpy array
            The T (thickness along s0) value (if None, it is internally calculated).

        Returns
        -------
        tuple
            the terms of the transfer matrix: (m11, m12, m21, m22).

        """

        if alpha is None or guigay_b is None or T is None:
            photon_out = self._calculatePhotonOut(photon_in)
            alpha = self._calculateAlphaGuigay(photon_in)
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
        Calculate the scattering matrix from the known transfer matrix (see equation XXX in Guigay and Sanchez del Rio).

        Parameters
        ----------
        transfer_matrix : tuple
            the terms of the transfer matrix: (m11, m12, m21, m22).

        Returns
        -------
        tuple
            the terms of the scattering matrix: (s11, s12, s21, s22).

        """
        m11, m12, m21, m22 = transfer_matrix
        s11 = m11 - m12 * m21 / m22
        s12 = m12 / m22
        s21 = -m21 / m22
        s22 = 1 / m22
        return s11, s12, s21, s22

    def calculateScatteringMatrix(self, photon_in, polarization=0, alpha=None, guigay_b=None, T=None):
        """
        Calculates the terms of the scattering matrix (see equation XXX in Guigay and Sanchez del Rio).
        Parameters
        ----------
        photon_in : instance of ComplexAmplitudePhoton
            Incoming photon or photon bunch.
        polarization : int
             0=sigma, 1=pi.
        is_thick : int
             Use thick crystal approximation: 0=No, 1=Yes.
        alpha : float or numpy array
            The alpha value (if None, it is internally calculated).
        guigay_b : float or numpy array
            The b (asymmetry factor) value (if None, it is internally calculated).
        T : float or numpy array
            The T (thickness along s0) value (if None, it is internally calculated).

        Returns
        -------
        tuple
            the terms of the scattering matrix: (s11, s12, s21, s22).

        """
        transfer_matrix = self.calculateTransferMatrix(photon_in, polarization=polarization,
                                                       alpha=alpha, guigay_b=guigay_b, T=T)
        return self.calculateScatteringMatrixFromTransferMatrix(transfer_matrix)


if __name__ == "__main__":
    a = CalculationStrategyMPMath()

    a0 = a.createVariable(numpy.array(0))
    print(type(a))
    print("cos(0)", a.cos(a0))

    api = a.createVariable(numpy.array(numpy.pi))
    print("cos(pi)", a.cos(api))

    api = a.createVariable(numpy.array([numpy.pi] * 10))
    print("cos(pi)", a.cos(api))
    # print("exp(pi)", a.exponentiate(api))

