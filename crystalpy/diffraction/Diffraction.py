"""
Defines a crystal diffraction experiment setup.
"""

from math import isnan
import numpy

from numpy import pi
import scipy.constants.codata

from crystalpy.diffraction.GeometryType import BraggDiffraction, BraggTransmission, LaueDiffraction, LaueTransmission
from crystalpy.diffraction.DiffractionExceptions import ReflectionImpossibleException, TransmissionImpossibleException, \
                                                                    StructureFactorF0isZeroException, StructureFactorFHisZeroException, \
                                                                    StructureFactorFHbarIsZeroException
from crystalpy.util.Photon import Photon
from crystalpy.util.PolarizedPhoton import PolarizedPhoton
from crystalpy.util.ComplexAmplitudePhotonBunch import ComplexAmplitudePhotonBunch

from crystalpy.util.PhotonBunch import PhotonBunch
from crystalpy.util.PolarizedPhotonBunch import PolarizedPhotonBunch


from crystalpy.diffraction.DiffractionResult import DiffractionResult
from crystalpy.diffraction.PerfectCrystalDiffraction import PerfectCrystalDiffraction
from crystalpy.polarization.CrystalPhasePlate import CrystalPhasePlate

from crystalpy.diffraction.DiffractionSetupSweeps import DiffractionSetupSweeps

class Diffraction(object):
    """This class contains methods to to make different diffraction calculations. It does not contain data."""

    #
    # main routines to calculate complex reflectivity and transmitivity
    #
    @classmethod
    def calculateDiffractedComplexAmplitudes(cls, diffraction_setup,
                                             incoming_photon,
                                             calculation_method=0,
                                             is_thick=0,
                                             use_transfer_matrix=0,
                                             calculation_strategy_flag=0,
                                             ):
        """Calculates the diffracted complex amplitude

        Parameters
        ----------
        diffraction_setup :

        incoming_photon :

        calculation_method : int, optional
             0: Zachariasen, 1: Guigay (Default value = 0)

        is_thick : int, optional
             0: No, 1: Yes (Default value = 0)

        use_transfer_matrix : int, optional
             0: No, 1: Yes (Default value = 0)

        calculation_strategy_flag : int, optional
            For computing exp, sin, cos:
            0: use mpmath, 1: use numpy, 2=use numpy truncated.

        Returns
        -------
        dict
            the complex amplitudes for sigma and pi polarizations:  {"S": float, "P": float}

        """
        # print(">>>> in calculateDiffractedComplexAmplitudes", use_transfer_matrix)

        # Get PerfectCrystal instance for the current photon.
        perfect_crystal = cls._perfectCrystalForPhoton(diffraction_setup,
                                                       incoming_photon,
                                                       calculation_strategy_flag=calculation_strategy_flag)

        # Calculate diffraction for current incoming photon.
        complex_amplitudes = perfect_crystal.calculateDiffraction(incoming_photon,
                                                                  calculation_method=calculation_method,
                                                                  is_thick=is_thick,
                                                                  use_transfer_matrix=use_transfer_matrix)

        return complex_amplitudes

    # using ComplexAmplitudePhoton
    @classmethod
    def calculateDiffractedComplexAmplitudePhoton(cls, diffraction_setup, photon,
                                                  calculation_method=0,
                                                  is_thick=0,
                                                  use_transfer_matrix=0,
                                                  calculation_strategy_flag=0,
                                                  ):
        """

        Parameters
        ----------
        diffraction_setup :

        photon :

        calculation_method : int, optional
             0: Zachariasen, 1: Guigay (Default value = 0)

        is_thick : int, optional
             0: No, 1: Yes (Default value = 0)

        use_transfer_matrix : int, optional
             0: No, 1: Yes (Default value = 0)

        calculation_strategy_flag : int, optional
            For computing exp, sin, cos:
            0: use mpmath, 1: use numpy, 2=use numpy truncated.

        Returns
        -------
        ComplexAmplitudePhoton instance

        """

        # print(">>>> in calculateDiffractedComplexAmplitudePhoton")

        # Get PerfectCrystal instance for the current photon.
        perfect_crystal = cls._perfectCrystalForPhoton(diffraction_setup,
                                                       photon,
                                                       calculation_strategy_flag=calculation_strategy_flag)

        coeffs = cls.calculateDiffractedComplexAmplitudes(diffraction_setup, photon,
                                                          calculation_method=calculation_method,
                                                          is_thick=is_thick,
                                                          use_transfer_matrix=use_transfer_matrix,
                                                          calculation_strategy_flag=calculation_strategy_flag)

        # Calculate outgoing Photon.
        outgoing_photon = perfect_crystal._calculatePhotonOut(photon)
        # apply reflectivities
        outgoing_photon.rescaleEsigma(coeffs["S"])
        outgoing_photon.rescaleEpi(coeffs["P"])

        return outgoing_photon

    @classmethod
    def calculateDiffractedComplexAmplitudePhotonBunch(cls,
                                                       diffraction_setup,
                                                       incoming_bunch,
                                                       calculation_method=0,
                                                       is_thick=0,
                                                       use_transfer_matrix=0,
                                                       calculation_strategy_flag=0):
        """Calculates the diffraction/transmission given by the setup.

        Parameters
        ----------
        diffraction_setup : DiffractionSetup instance
            The diffraction setup.

        incoming_bunch : ComplexAmplitudePhotonBeam instance

        calculation_method : int, optional
             0: Zachariasen, 1: Guigay (Default value = 0)

        is_thick : int, optional
             0: No, 1: Yes (Default value = 0)

        use_transfer_matrix : int, optional
             0: No, 1: Yes (Default value = 0)

        calculation_strategy_flag : int, optional
            For computing exp, sin, cos:
            0: use mpmath, 1: use numpy, 2=use numpy truncated.

        Returns
        -------
        ComplexAmplitudePhotonBeam instance
            Bunch made of diffracted/transmitted photons.

        """

        # print(">>>> in calculateDiffractedComplexAmplitudePhotonBunch", calculation_method)

        # Check that photon_bunch is indeed a PhotonBunch object.
        if not isinstance(incoming_bunch, ComplexAmplitudePhotonBunch):
            raise Exception("The incoming photon bunch must be a ComplexAmplitudePhotonBunch object!")

        vectorized_method = 2 # todo: delete options 0,1

        if vectorized_method == 2:
            return cls.calculateDiffractedComplexAmplitudePhoton(diffraction_setup,
                                                          incoming_bunch,
                                                          calculation_method=calculation_method,
                                                          is_thick=0,
                                                          use_transfer_matrix=0,
                                                          calculation_strategy_flag=calculation_strategy_flag)
        elif vectorized_method == 1:
            perfect_crystal = cls._perfectCrystalForPhotonBunch(diffraction_setup, incoming_bunch)
            coeffs = perfect_crystal.calculateDiffraction(incoming_bunch,
                                                          calculation_method=calculation_method,
                                                          is_thick=is_thick,
                                                          use_transfer_matrix=use_transfer_matrix,
                                                          calculation_strategy_flag=calculation_strategy_flag)
            outgoing_bunch = perfect_crystal._calculatePhotonOut(incoming_bunch)
            outgoing_bunch.rescaleEsigma(coeffs["S"])
            outgoing_bunch.rescaleEpi(coeffs["P"])
            # Return diffraction results.
            return outgoing_bunch
        else:
            # Create PhotonBunch instance.
            outgoing_bunch = ComplexAmplitudePhotonBunch([])

            perfect_crystal_bunch = cls._perfectCrystalForPhotonBunch(diffraction_setup,
                                                                      incoming_bunch,
                                                                      calculation_strategy_flag=calculation_strategy_flag)
            outgoing_bunch2 = perfect_crystal_bunch._calculatePhotonOut(incoming_bunch)

            for index, complex_amplitude_photon in enumerate(incoming_bunch):
                # Get PerfectCrystal instance for the current photon.
                perfect_crystal = cls._perfectCrystalForPhoton(diffraction_setup,
                                                               complex_amplitude_photon,
                                                               calculation_strategy_flag=calculation_strategy_flag)
                coeffs = perfect_crystal.calculateDiffraction(complex_amplitude_photon,
                                                              calculation_method=calculation_method,
                                                              is_thick=is_thick,
                                                              use_transfer_matrix=use_transfer_matrix,
                                                              calculation_strategy_flag=calculation_strategy_flag)

                # Calculate outgoing Photon.
                outgoing_complex_amplitude_photon = perfect_crystal._calculatePhotonOut(complex_amplitude_photon)
                # outgoing_complex_amplitude_photon = outgoing_bunch2.getPhotonIndex(index)
                # apply reflectivities
                outgoing_complex_amplitude_photon.rescaleEsigma(coeffs["S"])
                outgoing_complex_amplitude_photon.rescaleEpi(coeffs["P"])

                # Add result of current deviation.
                outgoing_bunch.addPhoton(outgoing_complex_amplitude_photon)

            # Return diffraction results.
            return outgoing_bunch

    #
    #
    #

    @classmethod
    def calculateDiffractedPolarizedPhoton(cls, diffraction_setup,
                                           incoming_polarized_photon,
                                           inclination_angle,
                                           calculation_method=0,
                                           is_thick=0,
                                           use_transfer_matrix=0,
                                           calculation_strategy_flag=0,
                                           ):
        """Calculates the diffraction/transmission given by the setup.

        Parameters
        ----------
        diffraction_setup : DiffractionSetup instance
            The diffraction setup.

        incoming_polarized_photon : ComplexAmplitudePhoton instance

        inclination_angle : float
            The inclination angle in rad.

        calculation_method : int, optional
             0: Zachariasen, 1: Guigay (Default value = 0)

        is_thick : int, optional
             0: No, 1: Yes (Default value = 0)

        use_transfer_matrix : int, optional
             0: No, 1: Yes (Default value = 0)

        calculation_strategy_flag : int, optional
            For computing exp, sin, cos:
            0: use mpmath, 1: use numpy, 2=use numpy truncated.

        Returns
        -------
        ComplexAmplitudePhoton instance
            Photon diffracted/transmitted.

        """

        # print(">>>> in calculateDiffractedPolarizedPhoton")

        # Retrieve the incoming Stokes vector.
        incoming_stokes_vector = incoming_polarized_photon.stokesVector()

        # Get PerfectCrystal instance for the current photon.
        perfect_crystal = cls._perfectCrystalForPhoton(diffraction_setup,
                                                       incoming_polarized_photon,
                                                       calculation_strategy_flag=calculation_strategy_flag)

        # Calculate diffraction for current incoming photon.
        complex_amplitudes = perfect_crystal.calculateDiffraction(incoming_polarized_photon,
                                                                  calculation_method=calculation_method,
                                                                  is_thick=is_thick,
                                                                  use_transfer_matrix=use_transfer_matrix,)

        # Calculate outgoing Photon.
        outgoing_photon = perfect_crystal._calculatePhotonOut(incoming_polarized_photon)

        # Calculate intensities and phases of the crystal  reflectivities or transmitivities
        intensity_pi = numpy.abs(complex_amplitudes["P"]) ** 2  # complex_amplitudes["P"].intensity()
        intensity_sigma = numpy.abs(complex_amplitudes["S"]) ** 2  # complex_amplitudes["S"].intensity()
        phase_pi = numpy.angle(numpy.array(complex_amplitudes["P"], dtype=complex))  # complex_amplitudes["P"].phase()
        phase_sigma = numpy.angle(
            numpy.array(complex_amplitudes["S"], dtype=complex))  # complex_amplitudes["S"].phase()

        # Get a CrystalPhasePlate instance which contains the Mueller matrix
        phase_plate = CrystalPhasePlate(  # incoming_stokes_vector=incoming_stokes_vector,
            intensity_sigma=intensity_sigma,
            phase_sigma=phase_sigma,
            intensity_pi=intensity_pi,
            phase_pi=phase_pi,
            inclination_angle=inclination_angle)

        # Use intensities and phases to calculate the Stokes vector for the outgoing photon.
        outgoing_stokes_vector = phase_plate.calculate_stokes_vector(incoming_stokes_vector)

        # Piece together the PolarizedPhoton object.
        outgoing_polarized_photon = PolarizedPhoton(energy_in_ev=outgoing_photon.energy(),
                                                    direction_vector=outgoing_photon.unitDirectionVector(),
                                                    stokes_vector=outgoing_stokes_vector)

        return outgoing_polarized_photon


    @classmethod
    def calculateDiffractedPolarizedPhotonBunch(cls,
                                                diffraction_setup,
                                                incoming_bunch,
                                                inclination_angle,
                                                calculation_method=0,
                                                is_thick=0,
                                                use_transfer_matrix=0,
                                                calculation_strategy_flag=0,
                                                ):
        """Calculates the diffraction/transmission bunch given by the crystal in the setup.

        Parameters
        ----------
        diffraction_setup : DiffractionSetup instance
            The diffraction setup.

        incoming_bunch : ComplexAmplitudePhotonBunch instance

        inclination_angle : float
            The inclination angle in rad.

        calculation_method : int, optional
             0: Zachariasen, 1: Guigay (Default value = 0)

        is_thick : int, optional
             0: No, 1: Yes (Default value = 0)

        use_transfer_matrix : int, optional
             0: No, 1: Yes (Default value = 0)

        calculation_strategy_flag : int, optional
            For computing exp, sin, cos:
            0: use mpmath, 1: use numpy, 2=use numpy truncated.

        Returns
        -------
        ComplexAmplitudePhotonBunch
            PhotonBunch object made up of diffracted/transmitted photons.

        """
        # print(">>>> in calculateDiffractedPolarizedPhotonBunch")

        # Create PhotonBunch instance.
        outgoing_bunch = PolarizedPhotonBunch([])

        # Check that photon_bunch is indeed a PhotonBunch object.
        if not isinstance(incoming_bunch, PolarizedPhotonBunch):
            raise Exception("The incoming photon bunch must be a PolarizedPhotonBunch object!")

        for index, polarized_photon in enumerate(incoming_bunch):
            outgoing_polarized_photon = cls.calculateDiffractedPolarizedPhoton(diffraction_setup, polarized_photon,
                                                                               inclination_angle,
                                                                               calculation_method=calculation_method,
                                                                               is_thick=is_thick,
                                                                               use_transfer_matrix=use_transfer_matrix,
                                                                               calculation_strategy_flag=calculation_strategy_flag)
            # Add result of current deviation.
            outgoing_bunch.addPhoton(outgoing_polarized_photon)

        # Return diffraction results.
        return outgoing_bunch

    @classmethod
    def calculateDiffraction(cls, diffraction_setup, calculation_method=0, is_thick=0, use_transfer_matrix=0,
                             calculation_strategy_flag=0):
        """Calculates the diffraction/transmission given by the setup.

        Parameters
        ----------
        diffraction_setup : DiffractionSetup instance

        calculation_method : int, optional
             0: Zachariasen, 1: Guigay (Default value = 0)

        is_thick : int, optional
             0: No, 1: Yes (Default value = 0)

        use_transfer_matrix : int, optional
             0: No, 1: Yes (Default value = 0)

        calculation_strategy_flag : int, optional
            For computing exp, sin, cos:
            0: use mpmath, 1: use numpy, 2=use numpy truncated.

        Returns
        -------
        DiffractionResult instance
            DiffractionResult with the input setup.

        """
        # print(">>>> in calculateDiffraction")
        if not isinstance(diffraction_setup, DiffractionSetupSweeps):
            raise Exception("Input object must be of type DiffractionSetupSweeps")

        # Create DiffractionResult instance.
        result = DiffractionResult(diffraction_setup, 0.0)

        for energy in diffraction_setup.energies():
            cls._calculateDiffractionForEnergy(diffraction_setup, energy, result,
                                               calculation_method=calculation_method,
                                               is_thick=is_thick,
                                               use_transfer_matrix=use_transfer_matrix,
                                               calculation_strategy_flag=calculation_strategy_flag)

        # Return diffraction results.
        return result

    #
    # checking routines
    #
    @classmethod
    def _checkSetup(cls, diffraction_setup, bragg_angle, F_0, F_H, F_H_bar):
        """Checks if a given diffraction setup is possible, i.e. if a given Diffraction/Transmission for the given asymmetry
        and Miller indices is possible. Raises an exception if impossible.

        Parameters
        ----------
        diffraction_setup :
            Diffraction setup.

        bragg_angle :
            Bragg angle.

        F_0 :
            Structure factor F_0.

        F_H :
            Structure factor F_H.

        F_H_bar :
            Structure factor F_H_bar.

        Raises
        ------
        Exception
            If the setup is not possible.

        """
        cls._checkSetupDiffraction(diffraction_setup, bragg_angle)
        cls._checkSetupStructureFactor(F_0, F_H, F_H_bar)

    @classmethod
    def _checkSetupStructureFactor(cls, F_0, F_H, F_H_bar):
        """Checks if the structure factor has reasonable values

        Parameters
        ----------
        F_0 :
            Structure factor F_0.

        F_H :
            Structure factor F_H.

        F_H_bar :
            Structure factor F_H_bar.

        Raises
        ------
        Exception
            If the structure factor values are not reasonable

        """

        # print(">>>> in _checkSetupStructureFactor")
        # Check structure factor F_0.
        if abs(F_0.real) < 1e-7 or isnan(F_0.real):
            raise StructureFactorF0isZeroException()

        # Check structure factor F_H.
        if abs(F_H.real) < 1e-7 or isnan(F_H.real) or abs(F_H.imag) < 1e-7 or isnan(F_H.imag):
            raise StructureFactorFHisZeroException()

        # Check structure factor F_H_bar.
        if abs(F_H_bar.real) < 1e-7 or isnan(F_H_bar.real) or abs(F_H_bar.imag) < 1e-7 or isnan(F_H_bar.imag):
            raise StructureFactorFHbarIsZeroException()

    @classmethod
    def _checkSetupDiffraction(cls, diffraction_setup, bragg_angle):
        """Checks if a given diffraction setup is possible, i.e. if a given Diffraction/Transmission for the given asymmetry
        and Miller indices is possible. Raises an exception if impossible.

        Parameters
        ----------
        diffraction_setup :
            Diffraction setup.

        bragg_angle :
            Bragg angle.

        Raises
        ------
        Exception
            If the this setup is not possible.

        """
        # print(">>>> in _checkSetupDiffraction")
        if (numpy.array(bragg_angle)).size > 1:
            asymmetry_angle = numpy.array(diffraction_setup.asymmetryAngle())
            # Check if the given geometry is a valid Bragg/Laue geometry.
            if diffraction_setup.geometryType() == BraggDiffraction() or diffraction_setup.geometryType() == BraggTransmission():
                if any(asymmetry_angle >= bragg_angle):
                    raise ReflectionImpossibleException()
            elif diffraction_setup.geometryType() == LaueDiffraction() or diffraction_setup.geometryType() == LaueTransmission():
                if any(asymmetry_angle <= bragg_angle):
                    raise TransmissionImpossibleException()
        else:
            asymmetry_angle = numpy.array(diffraction_setup.asymmetryAngle())
            # Check if the given geometry is a valid Bragg/Laue geometry.
            if diffraction_setup.geometryType() == BraggDiffraction() or diffraction_setup.geometryType() == BraggTransmission():
                if (asymmetry_angle >= bragg_angle):
                    raise ReflectionImpossibleException()
            elif diffraction_setup.geometryType() == LaueDiffraction() or diffraction_setup.geometryType() == LaueTransmission():
                if (asymmetry_angle <= bragg_angle):
                    raise TransmissionImpossibleException()




    # ##################################################################################################
    # FUNCTIONS ADAPTED TO WORK WITH A PHOTON OR PHOTON BUNCH
    # ##################################################################################################

    @classmethod
    def _perfectCrystalForEnergy(cls, diffraction_setup, energy,
                                 geometry_type=None,
                                 bragg_normal=None,
                                 surface_normal=None,
                                 # bragg_angle=None,
                                 # psi_0=None,
                                 # psi_H=None,
                                 # psi_H_bar=None,
                                 thickness=None,
                                 d_spacing=None,
                                 calculation_strategy_flag=0,
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

        calculation_strategy_flag : int, optional
            For computing exp, sin, cos:
            0: use mpmath, 1: use numpy, 2=use numpy truncated.

        Returns
        -------
        PerfectCrystalDiffraction instance

        """

        return PerfectCrystalDiffraction.initializeFromDiffractionSetupAndEnergy(diffraction_setup, energy,
                                                    geometry_type=geometry_type,
                                                    bragg_normal=bragg_normal,
                                                    surface_normal=surface_normal,
                                                    # bragg_angle=None,
                                                    # psi_0=None,
                                                    # psi_H=None,
                                                    # psi_H_bar=None,
                                                    thickness=thickness,
                                                    d_spacing=d_spacing,
                                                    calculation_strategy_flag=calculation_strategy_flag,
                                                    )

        # #
        # # energy-depending variables
        # #
        #
        # # Retrieve bragg angle.
        # bragg_angle = diffraction_setup.angleBragg(energy)
        #
        # # Check if given Bragg/Laue geometry and given miller indices are possible.
        # cls._checkSetupDiffraction(diffraction_setup, bragg_angle)
        #
        # psi_0, psi_H, psi_H_bar = diffraction_setup.psiAll(energy)
        #
        #
        # # Create PerfectCrystalDiffraction instance.
        # perfect_crystal = PerfectCrystalDiffraction(
        #     geometry_type   = geometry_type  if not geometry_type  is None else diffraction_setup.geometryType(),
        #     bragg_normal    = bragg_normal   if not bragg_normal   is None else diffraction_setup.vectorH(),
        #     surface_normal  = surface_normal if not surface_normal is None else diffraction_setup.vectorNormalSurface(),
        #     bragg_angle     = bragg_angle,
        #     psi_0           = psi_0,
        #     psi_H           = psi_H,
        #     psi_H_bar       = psi_H_bar,
        #     thickness       = thickness      if not thickness      is None else diffraction_setup.thickness(),
        #     d_spacing       = d_spacing      if not d_spacing      is None else diffraction_setup.dSpacing() * 1e-10,
        # )
        #
        # return perfect_crystal

    @classmethod
    def _perfectCrystalForPhoton(cls, diffraction_setup, polarized_photon, calculation_strategy_flag=0,):
        """

        Parameters
        ----------
        diffraction_setup :
            
        polarized_photon :

        calculation_strategy_flag : int, optional
            For computing exp, sin, cos:
            0: use mpmath, 1: use numpy, 2=use numpy truncated.

        Returns
        -------
        PerfectCrystalDiffraction instance

        """

        return cls._perfectCrystalForEnergy(diffraction_setup,
                                            polarized_photon.energy(),
                                            calculation_strategy_flag=calculation_strategy_flag)


    @classmethod
    def _perfectCrystalForPhotonBunch(cls, diffraction_setup, incoming_bunch, calculation_strategy_flag=0):
        """

        Parameters
        ----------
        diffraction_setup :
            
        incoming_bunch :

        calculation_strategy_flag : int, optional
            For computing exp, sin, cos:
            0: use mpmath, 1: use numpy, 2=use numpy truncated.

        Returns
        -------
        PerfectCrystalDiffraction instance

        """
        return cls._perfectCrystalForPhoton(diffraction_setup,
                                            incoming_bunch,
                                            calculation_strategy_flag=calculation_strategy_flag)


    # these methods use DiffractionSetupSweeps (for scans)
    # used in Diffraction.calculateDiffraction()
    @classmethod
    def _calculateDiffractionForEnergy(cls,
                                       diffraction_setup,
                                       energy,
                                       result,
                                       calculation_method=0,
                                       is_thick=0,
                                       use_transfer_matrix=0,
                                       calculation_strategy_flag=0):
        """Calculates the diffraction/transmission given by the setup.

        Parameters
        ----------
        diffraction_setup : DiffractionSetup instance
            The diffraction setup.

        energy : float

        result :  DiffractionResult instance, where results are added.

        calculation_method : int, optional
             0: Zachariasen, 1: Guigay (Default value = 0)

        is_thick : int, optional
             0: No, 1: Yes (Default value = 0)

        use_transfer_matrix : int, optional
             0: No, 1: Yes (Default value = 0)

        Returns
        -------
        DiffractionResult instance
            DiffractionResult with the calculation for the input energy added.

        """
        # Get PerfectCrystal instance for the current energy.
        if not isinstance(diffraction_setup, DiffractionSetupSweeps):
            raise Exception("Inmut must be of type: DiffractionSetupSweeps")

        perfect_crystal = cls._perfectCrystalForEnergy(diffraction_setup,
                                                       energy,
                                                       calculation_strategy_flag=calculation_strategy_flag)

        # For every deviation from Bragg angle ...
        for index, deviation in enumerate(diffraction_setup.angleDeviationGrid()):
            # Calculate deviated incoming photon.
            photon_direction = diffraction_setup.vectorIncomingPhotonDirection(energy, deviation)
            photon_in = Photon(energy, photon_direction)

            # Calculate diffraction for current incoming photon.
            result_deviation = perfect_crystal.calculateDiffraction(photon_in,
                                                                    calculation_method=calculation_method,
                                                                    is_thick=is_thick,
                                                                    use_transfer_matrix=use_transfer_matrix)

            # Calculate polarization difference between pi and sigma polarization.
            polarization_difference = result_deviation["P"] / result_deviation["S"]

            # Add result of current deviation.
            result.add(energy,
                       deviation,
                       result_deviation["S"],
                       result_deviation["P"],
                       polarization_difference)

        # Return diffraction results.
        return result


if __name__ == "__main__":
    if False:
        import numpy
        from crystalpy.diffraction.DiffractionSetupXraylib import DiffractionSetupXraylib
        from crystalpy.diffraction.GeometryType import BraggDiffraction, LaueDiffraction
        from crystalpy.diffraction.Diffraction import Diffraction
        from crystalpy.util.Vector import Vector


        diffraction_setup_r = DiffractionSetupXraylib(geometry_type=BraggDiffraction(),  # GeometryType object
                                               crystal_name="Si",  # string
                                               thickness=100e-6,  # meters
                                               miller_h=1,  # int
                                               miller_k=1,  # int
                                               miller_l=1,  # int
                                               asymmetry_angle=0,  # 10.0*numpy.pi/180.,            # radians
                                               azimuthal_angle=0.0)  # radians                            # int

        energy_setup = 8000.0
        bragg_angle = diffraction_setup_r.angleBragg(energy_setup)
        print("Bragg angle for E=%f eV is %f deg" % (energy_setup, bragg_angle * 180.0 / numpy.pi))

        deviation = 3e-6  # angle_deviation_min + ia * angle_step
        angle = deviation + bragg_angle

        # calculate the components of the unitary vector of the incident photon scan. Note that diffraction plane is YZ
        photon = Photon(energy_in_ev=energy_setup, direction_vector=Vector(0.0,
                                                                           numpy.cos(angle),
                                                                           - numpy.abs(numpy.sin(angle))))


        diffraction = Diffraction()
        coeffs_r = diffraction.calculateDiffractedComplexAmplitudes(diffraction_setup_r, photon)
        print(coeffs_r['S'], coeffs_r['P'])


        diffraction1 = Diffraction()
        coeffs_r = diffraction1.calculateDiffractedComplexAmplitudes(diffraction_setup_r, photon)
        print(coeffs_r['S'], coeffs_r['P'])

