"""
Represents Mueller calculation results.
"""
import numpy

# TODO inheritate from DiffractionResults?
class MuellerResult(object):
    """Constructor.

    Parameters
    ----------
    diffraction_result : DiffractionResult instance
        result of the diffraction.

    See Also
    --------
    crystalpy.diffraction.DiffractionResults.DiffractionResults

    """
    def __init__(self, diffraction_result):

        self.diffraction_result = diffraction_result
        self.diffraction_setup = diffraction_result.diffractionSetup()

        number_energies = len(self.energies())
        number_angles = len(self.angle_deviations())

        # Stokes parameters.
        self._s0 = numpy.zeros((number_energies,
                                number_angles))

        self._s1 = numpy.zeros((number_energies,
                                number_angles))

        self._s2 = numpy.zeros((number_energies,
                                number_angles))

        self._s3 = numpy.zeros((number_energies,
                                number_angles))

        # degree of circular polarization.
        self._circular_polarization_degree = numpy.zeros((number_energies,
                                                 number_angles))

    def energies(self):
        """Returns the energies used for these results.

        Returns
        -------
        numpy array
            The array with energies in eV

        """
        return self.diffraction_result.energies()

    def _energy_index(self, energy):
        """Returns the index of the entry in the energies list that is closest to the given energy.

        Parameters
        ----------
        energy :
            Energy to find index for.

        Returns
        -------
        type
            Energy index that corresponds to the energy.

        """
        energy_index = abs(self.energies()-energy).argmin()
        return energy_index

    def angle_deviations(self):
        """Returns the angle deviations used for these results.

        Returns
        -------
        numpy array
            Angle deviations used for these results.

        """
        return self.diffraction_result.angleDeviations()

    def _deviation_index(self, deviation):
        """Returns the index of the entry in the angle deviations list that is closest to the given deviation.

        Parameters
        ----------
        deviation :
            Deviation to find index for.

        Returns
        -------
        int
            Deviation index that corresponds to the deviation.

        """
        deviation_index = abs(self.angle_deviations()-deviation).argmin()
        return deviation_index

    def s0_by_energy(self, energy):
        """Returns the S0 Stokes parameter.

        Parameters
        ----------
        energy : float
            Energy corresponding to the returned S0.

        Returns
        -------
        Stokesvector instance.
            S0.

        """
        energy_index = self._energy_index(energy)
        return self._s0[energy_index, :]

    def s1_by_energy(self, energy):
        """Returns the S1 Stokes parameter.

        Parameters
        ----------
        energy : float
            Energy corresponding to the returned S1.

        Returns
        -------
        Stokesvector instance.
            S1.

        """
        energy_index = self._energy_index(energy)
        return self._s1[energy_index, :]

    def s2_by_energy(self, energy):
        """Returns the S2 Stokes parameter.

        Parameters
        ----------
        energy : float
            Energy corresponding to the returned S2.

        Returns
        -------
        Stokesvector instance.
            S2.

        """
        energy_index = self._energy_index(energy)
        return self._s2[energy_index, :]

    def s3_by_energy(self, energy):
        """Returns the S3 Stokes parameter.

        Parameters
        ----------
        energy : float
            Energy corresponding to the returned S3.

        Returns
        -------
        Stokesvector instance.
            S3.

        """
        energy_index = self._energy_index(energy)
        return self._s3[energy_index, :]

    def polarization_degree_by_energy(self, energy):
        """Returns the degree of circular polarization.

        Parameters
        ----------
        energy : float
            Energy corresponding to the returned circular polarization value.

        Returns
        -------
        float
            degree of circular polarization.

        """
        energy_index = self._energy_index(energy)
        return self._circular_polarization_degree[energy_index, :]

    def s0_by_deviation(self, deviation):
        """Returns the S0 Stokes parameter for a given deviation value.

        Parameters
        ----------
        deviation : float
            Deviation corresponding to the returned S0.

        Returns
        -------
        float
            S0.

        """
        deviation_index = self._deviation_index(deviation)
        return self._s0[deviation_index, :]

    def s1_by_deviation(self, deviation):
        """Returns the S1 Stokes parameter for a given deviation value.

        Parameters
        ----------
        deviation : float
            Deviation corresponding to the returned S1.

        Returns
        -------
        float
            S1.

        """
        deviation_index = self._deviation_index(deviation)
        return self._s1[deviation_index, :]

    def s2_by_deviation(self, deviation):
        """Returns the S2 Stokes parameter for a given deviation value.

        Parameters
        ----------
        deviation : float
            Deviation corresponding to the returned S2.

        Returns
        -------
        float
            S2.

        """
        deviation_index = self._deviation_index(deviation)
        return self._s2[deviation_index, :]

    def s3_by_deviation(self, deviation):
        """Returns the S0 Stokes parameter for a given deviation value.

        Parameters
        ----------
        deviation : float
            Deviation corresponding to the returned S3.

        Returns
        -------
        float
            S3.

        """
        deviation_index = self._deviation_index(deviation)
        return self._s3[deviation_index, :]

    def polarization_degree_by_deviation(self, deviation):
        """Returns the degree of circular polarization for a given deviation.

        Parameters
        ----------
        deviation :
            Deviation corresponding to the returned degree of circular polarization.

        Returns
        -------
        float
            degree of circular polarization.

        """
        deviation_index = self._deviation_index(deviation)
        return self._circular_polarization_degree[deviation_index, :]

    def add(self, energy, deviation, stokes_vector):
        """Adds a diffraction result for a given energy and deviation.

        Parameters
        ----------
        energy : float
            The energy in eV

        deviation : float
            The deviation angle in rad.

        stokes_vector : StokesVector instance
            The stokes vector

        """
        energy_index = self._energy_index(energy)
        deviation_index = self._deviation_index(deviation)

        self._s0[energy_index, deviation_index] = stokes_vector.s0
        self._s1[energy_index, deviation_index] = stokes_vector.s1
        self._s2[energy_index, deviation_index] = stokes_vector.s2
        self._s3[energy_index, deviation_index] = stokes_vector.s3
        self._circular_polarization_degree[energy_index, deviation_index] = stokes_vector.circularPolarizationDegree()
