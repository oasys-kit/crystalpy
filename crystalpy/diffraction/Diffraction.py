"""
Calculates a crystal diffraction.
Except for energy all units are in SI. Energy is in eV.
"""

from math import isnan

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

    @classmethod
    def _checkSetup(cls, diffraction_setup, bragg_angle, F_0, F_H, F_H_bar):
        """
        Checks if a given diffraction setup is possible, i.e. if a given Diffraction/Transmission for the given asymmetry
        and Miller indices is possible. Raises an exception if impossible.
        :param diffraction_setup: Diffraction setup.
        :param bragg_angle: Bragg angle.
        :param F_0: Structure factor F_0.
        :param F_H: Structure factor F_H.
        :param F_H_bar: Structure factor F_H_bar.
        """
        cls._checkSetupDiffraction(diffraction_setup, bragg_angle)
        cls._checkSetupStructureFactor(F_0, F_H, F_H_bar)

    @classmethod
    def _checkSetupStructureFactor(cls, F_0, F_H, F_H_bar):
        """
        Checks if the structure factor has reasonable values
        :param F_0: Structure factor F_0.
        :param F_H: Structure factor F_H.
        :param F_H_bar: Structure factor F_H_bar.
        """

        print(">>>> in _checkSetupStructureFactor")
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
        """
        Checks if a given diffraction setup is possible, i.e. if a given Diffraction/Transmission for the given asymmetry
        and Miller indices is possible. Raises an exception if impossible.
        :param diffraction_setup: Diffraction setup.
        :param bragg_angle: Bragg angle.
        """

        # Check if the given geometry is a valid Bragg/Laue geometry.
        if diffraction_setup.geometryType() == BraggDiffraction() or diffraction_setup.geometryType() == BraggTransmission():
            if diffraction_setup.asymmetryAngle() >= bragg_angle:
                raise ReflectionImpossibleException()
        elif diffraction_setup.geometryType() == LaueDiffraction() or diffraction_setup.geometryType() == LaueTransmission():
            if diffraction_setup.asymmetryAngle() <= bragg_angle:
                raise TransmissionImpossibleException()

    @classmethod
    def _perfectCrystalForEnergy(cls, diffraction_setup, energy):

        # Retrieve bragg angle.
        angle_bragg = diffraction_setup.angleBragg(energy)

        # Check if given Bragg/Laue geometry and given miller indices are possible.
        cls._checkSetupDiffraction(diffraction_setup, angle_bragg)


        # Retrieve lattice spacing d.
        d_spacing = diffraction_setup.dSpacing() * 1e-10

        # Calculate the Bragg normal B_H.
        normal_bragg = diffraction_setup.vectorH()

        # Calculate the surface normal n.
        normal_surface = diffraction_setup.vectorNormalSurface()


        psi_0, psi_H, psi_H_bar = diffraction_setup.psiAll(energy)
        # Create PerfectCrystalDiffraction instance.
        perfect_crystal = PerfectCrystalDiffraction(geometry_type=diffraction_setup.geometryType(),
                                                    bragg_normal=normal_bragg,
                                                    surface_normal=normal_surface,
                                                    bragg_angle=angle_bragg,
                                                    psi_0=psi_0,
                                                    psi_H=psi_H,
                                                    psi_H_bar=psi_H_bar,
                                                    thickness=diffraction_setup.thickness(),
                                                    d_spacing=d_spacing)

        return perfect_crystal


    # calculate complex reflectivity and transmitivity
    @classmethod
    def calculateDiffractedComplexAmplitudes(cls, diffraction_setup,
                                             incoming_photon,
                                             calculation_method=0,
                                             is_thick=0,
                                             use_transfer_matrix=0,
                                             ):

        # Get PerfectCrystal instance for the current photon.
        perfect_crystal = cls._perfectCrystalForPhoton(diffraction_setup, incoming_photon)

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
                                                  ):

        # Get PerfectCrystal instance for the current photon.
        perfect_crystal = cls._perfectCrystalForPhoton(diffraction_setup, photon)

        coeffs = cls.calculateDiffractedComplexAmplitudes(diffraction_setup, photon,
                                                          calculation_method=calculation_method,
                                                          is_thick=is_thick,
                                                          use_transfer_matrix=use_transfer_matrix)

        # Calculate outgoing Photon.
        outgoing_photon = perfect_crystal._calculatePhotonOut(photon)
        # apply reflectivities
        outgoing_photon.rescaleEsigma(coeffs["S"])
        outgoing_photon.rescaleEpi(coeffs["P"])

        return outgoing_photon


    # ##################################################################################################
    # FUNCTIONS ADAPTED TO WORK WITH A PHOTON OR PHOTON BUNCH
    # ##################################################################################################

    @classmethod
    def _perfectCrystalForPhoton(cls, diffraction_setup, polarized_photon):

        energy = polarized_photon.energy()

        # Retrieve bragg angle.
        angle_bragg = diffraction_setup.angleBragg(energy)


        # Check if given Bragg/Laue geometry and given miller indices are possible.
        cls._checkSetupDiffraction(diffraction_setup, angle_bragg)


        # Retrieve lattice spacing d.
        d_spacing = diffraction_setup.dSpacing() * 1e-10

        # Calculate the Bragg normal B_H.
        normal_bragg = diffraction_setup.vectorH()

        # Calculate the surface normal n.
        normal_surface = diffraction_setup.vectorNormalSurface()

        psi_0, psi_H, psi_H_bar = diffraction_setup.psiAll(energy)

        # Create PerfectCrystalDiffraction instance.
        perfect_crystal = PerfectCrystalDiffraction(geometry_type=diffraction_setup.geometryType(),
                                                    bragg_normal=normal_bragg,
                                                    surface_normal=normal_surface,
                                                    bragg_angle=angle_bragg,
                                                    psi_0=psi_0,
                                                    psi_H=psi_H,
                                                    psi_H_bar=psi_H_bar,
                                                    thickness=diffraction_setup.thickness(),
                                                    d_spacing=d_spacing)

        return perfect_crystal

    @classmethod
    def _perfectCrystalForPhotonBunch(cls, diffraction_setup, incoming_bunch):

        energies = incoming_bunch.energies()

        # Retrieve bragg angle.
        angle_bragg = diffraction_setup.angleBragg(energies)

        # Check if given Bragg/Laue geometry and given miller indices are possible.
        cls._checkSetupDiffraction(diffraction_setup, angle_bragg[0])

        # Retrieve lattice spacing d.
        d_spacing = diffraction_setup.dSpacing() * 1e-10

        # Calculate the Bragg normal B_H.
        normal_bragg = diffraction_setup.vectorH()

        # Calculate the surface normal n.
        normal_surface = diffraction_setup.vectorNormalSurface()

        psi_0, psi_H, psi_H_bar = diffraction_setup.psiAll(energies)

        # Create PerfectCrystalDiffraction instance.
        perfect_crystal = PerfectCrystalDiffraction(geometry_type=diffraction_setup.geometryType(),
                                                    bragg_normal=normal_bragg,
                                                    surface_normal=normal_surface,
                                                    bragg_angle=angle_bragg,
                                                    psi_0=psi_0,
                                                    psi_H=psi_H,
                                                    psi_H_bar=psi_H_bar,
                                                    thickness=diffraction_setup.thickness(),
                                                    d_spacing=d_spacing)

        return perfect_crystal

    @classmethod
    def calculateDiffractedComplexAmplitudePhotonBunch(cls, diffraction_setup, incoming_bunch,
                                                       calculation_method=0,
                                                       is_thick=0,
                                                       use_transfer_matrix=0):
        """
        Calculates the diffraction/transmission given by the setup.
        :param diffraction_setup: The diffraction setup.
        :return: PhotonBunch object made up of diffracted/transmitted photons.
        """
        # Create PhotonBunch instance.
        outgoing_bunch = ComplexAmplitudePhotonBunch([])

        # Retrieve the photon bunch from the diffraction setup.
        # incoming_bunch = diffraction_setup.incomingPhotons()

        # Check that photon_bunch is indeed a PhotonBunch object.
        if not isinstance(incoming_bunch, ComplexAmplitudePhotonBunch):
            raise Exception("The incoming photon bunch must be a ComplexAmplitudePhotonBunch object!")

        vectorized_method = 1

        if vectorized_method:
            perfect_crystal = cls._perfectCrystalForPhotonBunch(diffraction_setup, incoming_bunch)
            coeffs = perfect_crystal.calculateDiffraction(incoming_bunch,
                                                          calculation_method=calculation_method,
                                                          is_thick=is_thick,
                                                          use_transfer_matrix=use_transfer_matrix)
            outgoing_bunch = perfect_crystal._calculatePhotonOut(incoming_bunch)
            outgoing_bunch.rescaleEsigma(coeffs["S"])
            outgoing_bunch.rescaleEpi(coeffs["P"])

        else:
            perfect_crystal_bunch = cls._perfectCrystalForPhotonBunch(diffraction_setup, incoming_bunch)
            outgoing_bunch2 = perfect_crystal_bunch._calculatePhotonOut(incoming_bunch)

            for index, complex_amplitude_photon in enumerate(incoming_bunch):

                # Get PerfectCrystal instance for the current photon.
                perfect_crystal = cls._perfectCrystalForPhoton(diffraction_setup, complex_amplitude_photon)
                coeffs = perfect_crystal.calculateDiffraction(complex_amplitude_photon,
                                                              calculation_method=calculation_method,
                                                              is_thick=is_thick,
                                                              use_transfer_matrix=use_transfer_matrix)

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


    @classmethod
    def calculateDiffractedPolarizedPhoton(cls, diffraction_setup, incoming_polarized_photon, inclination_angle,
                                           calculation_method=0,
                                           is_thick=0,
                                           use_transfer_matrix=0):
        """
        Calculates the diffraction/transmission given by the setup.
        :param diffraction_setup: The diffraction setup.
        :return: PhotonBunch object made up of diffracted/transmitted photons.
        """
        # Retrieve the incoming Stokes vector.
        incoming_stokes_vector = incoming_polarized_photon.stokesVector()

        # Get PerfectCrystal instance for the current photon.
        perfect_crystal = cls._perfectCrystalForPhoton(diffraction_setup, incoming_polarized_photon)

        # Calculate diffraction for current incoming photon.
        complex_amplitudes = perfect_crystal.calculateDiffraction(incoming_polarized_photon,
                                                                  calculation_method=calculation_method,
                                                                  is_thick=is_thick,
                                                                  use_transfer_matrix=use_transfer_matrix)

        # Calculate outgoing Photon.
        outgoing_photon = perfect_crystal._calculatePhotonOut(incoming_polarized_photon)

        # Calculate intensities and phases of the crystal  reflectivities or transmitivities
        intensity_pi = complex_amplitudes["P"].intensity()
        intensity_sigma = complex_amplitudes["S"].intensity()
        phase_pi = complex_amplitudes["P"].phase()
        phase_sigma = complex_amplitudes["S"].phase()

        # Get a CrystalPhasePlate instance which contains the Mueller matrix
        phase_plate = CrystalPhasePlate( #incoming_stokes_vector=incoming_stokes_vector,
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
    def calculateDiffractedPolarizedPhotonBunch(cls, diffraction_setup, incoming_bunch, inclination_angle,
                                                calculation_method=0,
                                                is_thick=0,
                                                use_transfer_matrix=0):
        """
        Calculates the diffraction/transmission given by the setup.
        :param diffraction_setup: The diffraction setup.
        :return: PhotonBunch object made up of diffracted/transmitted photons.
        """
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
                                                                               use_transfer_matrix=use_transfer_matrix)
            # Add result of current deviation.
            outgoing_bunch.addPhoton(outgoing_polarized_photon)

        # Return diffraction results.
        return outgoing_bunch

    # ##################################################################################################
    # these methods use DiffractionSetupSweeps (for scans)
    # ##################################################################################################
    @classmethod
    def _calculateDiffractionForEnergy(cls, diffraction_setup, energy, result,
                                       calculation_method=0,
                                       is_thick=0,
                                       use_transfer_matrix=0):
        """
        Calculates the diffraction/transmission given by the setup.
        :param diffraction_setup: The diffraction setup.
        :return: DiffractionResult representing this setup.
        """
        # Get PerfectCrystal instance for the current energy.
        if not isinstance(diffraction_setup, DiffractionSetupSweeps):
            raise Exception("Inmut must be of type: DiffractionSetupSweeps")

        perfect_crystal = cls._perfectCrystalForEnergy(diffraction_setup, energy)

        # For every deviation from Bragg angle ...
        for index, deviation in enumerate(diffraction_setup.angleDeviationGrid()):
            # Calculate deviated incoming photon.
            photon_direction = diffraction_setup.incomingPhotonDirection(energy, deviation)
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

    @classmethod
    def calculateDiffraction(cls, diffraction_setup, calculation_method=0, is_thick=0, use_transfer_matrix=0):
        """
        Calculates the diffraction/transmission given by the setup.
        :param diffraction_setup: The diffraction setup.
        :calculation_method: 0=Zachariasen, 1=Guigay
        :return: DiffractionResult representing this setup.
        """

        if not isinstance(diffraction_setup, DiffractionSetupSweeps):
            raise Exception("Input object must be of type DiffractionSetupSweeps")

        # Create DiffractionResult instance.
        result = DiffractionResult(diffraction_setup, 0.0)

        for energy in diffraction_setup.energies():
            cls._calculateDiffractionForEnergy(diffraction_setup, energy, result,
                                               calculation_method=calculation_method,
                                               is_thick=is_thick,
                                               use_transfer_matrix=use_transfer_matrix)

        # Return diffraction results.
        return result

if __name__ == "__main__":
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

