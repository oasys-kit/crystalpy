"""
Object to hold the diffraction results.

    self._intensities = numpy.array((number_energies,number_angles,index_polarization))
    self._phases = numpy.array((number_energies,number_angles,index_polarization))

    index_polarization:
        INDEX_POLARIZATION_S = 0
        INDEX_POLARIZATION_P = 1
        INDEX_DIFFERENCE_PS = 2

    Note that INDEX_DIFFERENCE_PS=2 means for:
        self._intensities: the RATIO of p-intensity / s-intensity
        self._phases: the DIFFERENCE of p-intensity - s-intensity

"""
import numpy


class DiffractionResult(object):
    """
    Constructor.

    Parameters
    ----------

    diffraction_setup : DiffractionSetup instance
        Setup used for these results.

    bragg_angle: float
        Bragg angle in rad of the setup.
    """

    INDEX_POLARIZATION_S = 0
    INDEX_POLARIZATION_P = 1
    INDEX_DIFFERENCE_PS = 2

    def __init__(self, diffraction_setup, bragg_angle):

        self._diffraction_setup = diffraction_setup.clone()
        self._bragg_angle = bragg_angle

        number_energies = len(self.energies())
        number_angles = len(self.angleDeviations())
        number_polarizations = 3

        self._intensities = numpy.zeros((number_energies,
                                         number_angles,
                                         number_polarizations))

        self._phases = numpy.zeros((number_energies,
                                    number_angles,
                                    number_polarizations))

    def diffractionSetup(self):
        """Returns the diffraction setup used for the calculation of these results.
        :return:

        Returns
        -------
        DiffractionSetup instance
            The diffraction setup used for the calculation of these results (referenced, not copied).
        """
        return self._diffraction_setup

    def braggAngle(self):
        """Returns Bragg angle used for these results.

        Returns
        -------
        float
            Bragg angle used for these results.

        """
        return self._bragg_angle

    def energies(self):
        """Returns the energies used for these results.

        Returns
        -------
        numpy array
            Energies used for these results.

        """
        return self._diffraction_setup.energies()

    def _energyIndexByEnergy(self, energy):
        """Returns the index of the entry in the energies list that is closest to the given energy.

        Parameters
        ----------
        energy : float
            Energy to find index for.

        Returns
        -------
        int
            Energy index that corresponds to the energy.

        """
        energy_index = abs(self.energies()-energy).argmin()
        return energy_index

    def angleDeviations(self):
        """Returns the angle deviations used for these results.

        Returns
        -------
        numpy array
            Angle deviations used for these results.

        """
        return self._diffraction_setup.angleDeviationGrid()

    def _deviationIndexByDeviation(self, deviation):
        """Returns the index of the entry in the angle deviations list that is closest to the given deviation.

        Parameters
        ----------
        deviation :
            Deviation corresponding to the index we are looking for.

        Returns
        -------
        int
            index that corresponds to the deviation.

        """
        deviation_index = abs(self.angleDeviations()-deviation).argmin()
        return deviation_index

    def angles(self):
        """Returns the angles used for calculation of these results.

        Returns
        -------
        numpy array
            The angles used for the calculation of these results.

        """
        return [self.braggAngle() + dev for dev in self.angleDeviations()]

    def sIntensityByEnergy(self, energy):
        """Returns the intensity of the S polarization.

        Parameters
        ----------
        energy : float
            Energy for the intensity we are loking for.

        Returns
        -------
        float
            Intensity of the S polarization.

        """
        energy_index = self._energyIndexByEnergy(energy)
        return self._intensities[energy_index, :, self.INDEX_POLARIZATION_S]

    def sPhaseByEnergy(self, energy, deg=False):
        """Returns the phase of the S polarization.

        Parameters
        ----------
        energy : float
            Energy for the phase we are loking for.

        deg : boolean, optional
            if True the phase is converted into degrees. (Default deg = False)

        Returns
        -------
        float
            Phase of the S polarization.

        """
        energy_index = self._energyIndexByEnergy(energy)
        if deg:
            return self._phases[energy_index, :, self.INDEX_POLARIZATION_S] * 180 / numpy.pi
        else:
            return self._phases[energy_index, :, self.INDEX_POLARIZATION_S]

    def pIntensityByEnergy(self, energy):
        """Returns the intensity of the P polarization.

        Parameters
        ----------
        energy : float
            Energy for the intensity we are loking for.

        Returns
        -------
        float
            Intensity of the S polarization.

        """
        energy_index = self._energyIndexByEnergy(energy)
        return self._intensities[energy_index, :, self.INDEX_POLARIZATION_P]

    def pPhaseByEnergy(self, energy, deg=False):
        """Returns the phase of the P polarization.

        Parameters
        ----------
        energy : float
            Energy for the phase we are loking for.

        deg : boolean, optional
            if True the phase is converted into degrees. (Default deg = False)

        Returns
        -------
        float
            Phase of the S polarization.

        """
        energy_index = self._energyIndexByEnergy(energy)
        if deg:
            return self._phases[energy_index, :, self.INDEX_POLARIZATION_P] * 180 / numpy.pi
        else:
            return self._phases[energy_index, :, self.INDEX_POLARIZATION_P]

    def differenceIntensityByEnergy(self, energy):
        """Returns the ratio of the intensity between P over S polarizations.

        Parameters
        ----------
        energy : float
            Energy to return intensity for.

        Returns
        -------
        float
            Ratio of the intensities P / S.

        """
        energy_index = self._energyIndexByEnergy(energy)
        return self._intensities[energy_index, :, self.INDEX_DIFFERENCE_PS]

    def differencePhaseByEnergy(self, energy, deg=False):
        """Returns the difference of phase between P and S polarizations.

        Parameters
        ----------
        energy : float
            Energy to return phase for.

        deg : boolean, optional
            if True the phase is converted into degrees. (Default deg = False)

        Returns
        -------
        float
            Phase of the difference between S and P polarization.

        """
        energy_index = self._energyIndexByEnergy(energy)
        if deg:
            return self._phases[energy_index, :, self.INDEX_DIFFERENCE_PS] * 180 / numpy.pi
        else:
            return self._phases[energy_index, :, self.INDEX_DIFFERENCE_PS]

    def sIntensityByDeviation(self, deviation):
        """Returns the intensity of the S polarization.

        Parameters
        ----------
        deviation : float
            Deviation to return intensity for.

        Returns
        -------
        float
            Intensity of the S polarization.

        """
        deviation_index = self._deviationIndexByDeviation(deviation)
        return self._intensities[:, deviation_index, self.INDEX_POLARIZATION_S]

    def sPhaseByDeviation(self, deviation, deg=False):
        """Returns the phase of the S polarization.

        Parameters
        ----------
        deviation : float
            Deviation to return phase for.

        deg : boolean, optional
            if True the phase is converted into degrees. (Default value = False)

        Returns
        -------
        float
            Phase of the S polarization.

        """
        deviation_index = self._deviationIndexByDeviation(deviation)
        if deg:
            return self._phases[deviation_index, :, self.INDEX_POLARIZATION_S] * 180 / numpy.pi
        else:
            return self._phases[:, deviation_index, self.INDEX_POLARIZATION_S]

    def pIntensityByDeviation(self, deviation):
        """Returns the intensity of the P polarization.

        Parameters
        ----------
        deviation : float
            Deviation to return intensity for.

        Returns
        -------
        float
            Intensity of the P polarization.

        """
        deviation_index = self._deviationIndexByDeviation(deviation)
        return self._intensities[:, deviation_index, self.INDEX_POLARIZATION_P]

    def pPhaseByDeviation(self, deviation, deg=False):
        """Returns the phase of the P polarization.

        Parameters
        ----------
        deviation : float
            Deviation to return phase for.

        deg : boolean, optional
            if True the phase is converted into degrees. (Default value = False)

        Returns
        -------
        float
            Phase of the P polarization.

        """
        deviation_index = self._deviationIndexByDeviation(deviation)
        if deg:
            return self._phases[deviation_index, :, self.INDEX_POLARIZATION_P] * 180 / numpy.pi
        else:
            return self._phases[:, deviation_index, self.INDEX_POLARIZATION_P]

    def differenceIntensityByDeviation(self, deviation):
        """Returns the intensities ratio for the two polarizations P/S .

        Parameters
        ----------
        deviation : float
            Deviation to return intensity for.

        Returns
        -------
        float
            Intensity ratio P / S for the two polarizations.

        """
        deviation_index = self._deviationIndexByDeviation(deviation)
        return self._intensities[:, deviation_index, self.INDEX_DIFFERENCE_PS]

    def differencePhaseByDeviation(self, deviation, deg=False):
        """Returns the difference of the phase between P and S polarizations.

        Parameters
        ----------
        deviation : float
            Deviation to return phase for.

        deg : boolean, optional
            if True the phase is converted into degrees. (Default value = False)

        Returns
        -------
        float
            Different of phase between P and S polarization.

        """
        deviation_index = self._deviationIndexByDeviation(deviation)
        if deg:
            return self._phases[deviation_index, :, self.INDEX_DIFFERENCE_PS] * 180 / numpy.pi
        else:
            return self._phases[:, deviation_index, self.INDEX_DIFFERENCE_PS]

    def add(self, energy, deviation, s_complex_amplitude, p_complex_amplitude, difference_complex_amplitude):
        """Adds a new result for a given energy and deviation.

        Parameters
        ----------
        energy : float
            The photon energy in eV.
            
        deviation : float
            The deviation angle in rad.
            
        s_complex_amplitude : float
            The complex amplitude for sigma polarization.
            
        p_complex_amplitude : float
            The complex amplitude for pi polarization.
            
        difference_complex_amplitude : float
            The ratio between the intensities of pi over sigma polarizations.

        """
        energy_index = self._energyIndexByEnergy(energy)
        deviation_index = self._deviationIndexByDeviation(deviation)

        self._intensities[energy_index, deviation_index, self.INDEX_POLARIZATION_S] = numpy.abs(s_complex_amplitude)**2
        self._intensities[energy_index, deviation_index, self.INDEX_POLARIZATION_P] = numpy.abs(p_complex_amplitude)**2
        self._intensities[energy_index, deviation_index, self.INDEX_DIFFERENCE_PS] = numpy.abs(difference_complex_amplitude)**2

        self._phases[energy_index, deviation_index, self.INDEX_POLARIZATION_S] = numpy.angle(numpy.array(s_complex_amplitude, dtype=complex))
        self._phases[energy_index, deviation_index, self.INDEX_POLARIZATION_P] = numpy.angle(numpy.array(p_complex_amplitude, dtype=complex))
        self._phases[energy_index, deviation_index, self.INDEX_DIFFERENCE_PS] = numpy.angle(numpy.array(difference_complex_amplitude, dtype=complex))

