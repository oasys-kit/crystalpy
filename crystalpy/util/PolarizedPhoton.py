"""
This object contains a list of PolarizedPhoton objects, characterized by energy, direction vector and Stokes vector.
This object is used as input to and output from the passive crystal widget.
"""

from crystalpy.util.Photon import Photon
from crystalpy.util.StokesVector import StokesVector
from crystalpy.polarization.MuellerMatrix import MuellerMatrix

class PolarizedPhoton(Photon):
    """A Photon object with a specified polarization state described by a Stokes vector."""
    def __init__(self, energy_in_ev, direction_vector, stokes_vector):
        """Constructor.

        Parameters
        ----------
        energy_in_ev : float
            Photon energy in eV.

        direction_vector : Vector instance
            The direction of the photon (no need to be normalized).

        stokes_vector : StokesVector instance
            Stokes vector describing the polarization state.

        Returns
        -------
        PolarizedPhoton instance


        """

        self._stokes_vector = stokes_vector
        super(PolarizedPhoton, self).__init__(energy_in_ev, direction_vector)

    def duplicate(self):
        """Duplicates a stokes photon.

        Returns
        -------
        Photon instance
            New StokesPhoton instance with identical photon.

        """
        return PolarizedPhoton(self._energy_in_ev,
                               self._unit_direction_vector.duplicate(),
                               self._stokes_vector.duplicate())


    def stokesVector(self):
        """Returns the Stokes vector."""
        return self._stokes_vector

    def setStokesVector(self,stokes_vector):
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
        mueller_matrix :
             (Default value = MuellerMatrix())

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
