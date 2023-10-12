"""
This object represents a polarized photon, characterized by energy, direction vector and Stokes parameters.

"""

from crystalpy.util.Photon import Photon
from crystalpy.util.StokesVector import StokesVector
from crystalpy.polarization.MuellerMatrix import MuellerMatrix

class PolarizedPhoton(Photon):
    """Constructor.

    Parameters
    ----------
    energy_in_ev : float
        Photon energy in eV.

    direction_vector : Vector instance
        The direction of the photon (no need to be normalized).

    stokes_vector : StokesVector instance
        Stokes vector describing the polarization state.

    """
    def __init__(self, energy_in_ev, direction_vector, stokes_vector):
        super(PolarizedPhoton, self).__init__(energy_in_ev, direction_vector)
        self._stokes_vector = stokes_vector

    # def duplicate(self):
    #     """Duplicates a stokes photon.
    #
    #     Returns
    #     -------
    #     PolarizedPhoton instance
    #         New PolarizedPhoton instance with identical photon.
    #
    #     """
    #     return PolarizedPhoton(self._energy_in_ev,
    #                            self._unit_direction_vector.duplicate(),
    #                            self._stokes_vector.duplicate())


    def stokesVector(self):
        """Returns the Stokes vector.

        Returns
        -------
        StokesVector instance
            The Stokes vector (referenced, not copied).
        """
        return self._stokes_vector

    def setStokesVector(self, stokes_vector):
        """Sets the stokes vector

        Parameters
        ----------
        stokes_vector : StokesVector instance

        """
        self._stokes_vector = stokes_vector

    def applyMuellerMatrix(self,mueller_matrix=MuellerMatrix()):
        """Modify the stokes vector by a Muller matrix.

        Parameters
        ----------
        mueller_matrix : MuellerMatrix instance, optional
            The Mueller matrix

        See Also
        --------
        crystalpy.polarization.MuellerMatrix.MuellerMatrix

        """
        s_in = self.stokesVector()
        s_out = mueller_matrix.calculate_stokes_vector( s_in )
        self.setStokesVector(s_out)


    def circularPolarizationDegree(self):
        """Returns the degree of circular polarization.

        Returns
        -------
        float
            The polarization degree.

        """
        return self._stokes_vector.circularPolarizationDegree()

    def __eq__(self, candidate):
        if ((self.energy() == candidate.energy() and
                self.unitDirectionVector() == candidate.unitDirectionVector()) and
                self.stokesVector() == candidate.stokesVector()):
            return True

        return False
