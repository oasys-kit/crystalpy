"""
Represents Mueller diffraction setup.
"""

from crystalpy.polarization.CrystalPhasePlate import CrystalPhasePlate
from crystalpy.polarization.MuellerResult import MuellerResult


class MuellerDiffraction(object):
    """Constructor.

    Parameters
    ----------
    diffraction_result : DiffractionResult instance
        result of the diffraction.

    incoming_stokes_vector : StokesVector instance
        the Stokes vectoR.

    inclination_angle : float, optional
        the inclination angle in rad.

    See Also
    --------
    crystalpy.diffraction.DiffractionResults.DiffractionResults

    """

    def __init__(self, diffraction_result, incoming_stokes_vector, inclination_angle=0.0):

        self._diffraction_result = diffraction_result  # DiffractionResult object.
        self._incoming_stokes_vector = incoming_stokes_vector  # StokesVector object.
        self._inclination_angle = inclination_angle  # radians.

    def _intensity_sigma(self, energy, index):
        """Retrieves the (forward-) diffracted intensity for the sigma polarization.

        Parameters
        ----------
        energy :
            energy for which the (forward-) diffracted intensity is calculated.
        index :
            

        Returns
        -------
        float
            intensity of the sigma polarization for the given index.

        """
        return self._diffraction_result.sIntensityByEnergy(energy)[index]

    def _phase_sigma(self, energy, index):
        """Retrieves the phase for the sigma polarization.

        Parameters
        ----------
        energy : float
            energy for which the phase is calculated.

        index : int
            the index of the wanted result.

        Returns
        -------
        float
            intensity of the sigma polarization.

        """
        return self._diffraction_result.sPhaseByEnergy(energy, deg=False)[index]

    def _intensity_pi(self, energy, index):
        """Retrieves the (forward-) diffracted intensity for the pi polarization.

        Parameters
        ----------
        energy :
            energy for which the (forward-) diffracted intensity is calculated.

        index : int
            the index of the wanted result.

        Returns
        -------
        float
            intensity of the pi polarization.

        """
        return self._diffraction_result.pIntensityByEnergy(energy)[index]

    def _phase_pi(self, energy, index):
        """Retrieves the phase for the pi polarization.

        Parameters
        ----------
        energy :
            energy for which the phase is calculated.

        index : int
            the index of the wanted result.

        Returns
        -------
        float
            intensity of the pi polarization.

        """
        return self._diffraction_result.pPhaseByEnergy(energy, deg=False)[index]

    # TODO: rename to add_photon??
    def _calculate_stokes_for_energy(self, energy, mueller_result):
        """Calculates the outgoing Stokes vectors (deviation) for a certain energy.
        :return: StokesVector objects (deviations).

        Parameters
        ----------
        energy : float
            The energy in eV
            
        mueller_result : MuellerResult instance
            The object where the results are added
            

        Returns
        -------

        """
        for index, deviation in enumerate(self._diffraction_result.angleDeviations()):

            intensity_sigma = self._intensity_sigma(energy, index)
            phase_sigma = self._phase_sigma(energy, index)
            intensity_pi = self._intensity_pi(energy, index)
            phase_pi = self._phase_pi(energy, index)

            crystal_phase_plate = CrystalPhasePlate( # self._incoming_stokes_vector,
                                                    intensity_sigma, phase_sigma,
                                                    intensity_pi, phase_pi,
                                                    self._inclination_angle)

            mueller_result.add(energy, deviation, crystal_phase_plate.calculate_stokes_vector(self._incoming_stokes_vector) )

    def calculate_stokes(self):
        """Calculates the outgoing Stokes vectors (deviation).

        Parameters
        ----------

        Returns
        -------
        MuellerResult instance


        """
        # Create an instance of the MuellerResult class.
        mueller_result = MuellerResult(self._diffraction_result)

        # Iterate over the different energies.
        for energy in self._diffraction_result.energies():
            self._calculate_stokes_for_energy(energy, mueller_result)

        return mueller_result
