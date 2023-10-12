"""
A stack of polarized photons characterized by photon energy, direction vector and Stokes vector.
"""
import numpy

from crystalpy.util.PhotonBunch import PhotonBunch, PhotonBunchDecorator
from crystalpy.util.PolarizedPhoton import PolarizedPhoton
from crystalpy.util.Vector import Vector
from crystalpy.util.StokesVector import StokesVector

#todo: delete
class PolarizedPhotonBunchOld(PhotonBunch):
    """The PolarizadPhotonBunch is is a collection of PolarizedPhoton objects, making up the polarized photon beam.

    Constructor.

    Parameters
    ----------
    polarized_photons : list, optional
        List of PolarizedPhoton instances.

    """


    def __init__(self, polarized_photons=None):
        if polarized_photons == None:
            self.polarized_photon_bunch = []
        else:
            self.polarized_photon_bunch = polarized_photons


    def toDictionary(self):
        """Created a dictionary containing information about the bunch.

        Returns
        -------
        dict
            Information in tags: "number of photons", "energies", "deviations", "vx", "vy", "vz", "s0", "s1", "s2", "s3" and "polarization degree".

        """

        array_dict = PhotonBunch.toDictionary(self)

        stokes = numpy.zeros([4, len(self)])
        polarization_degrees = numpy.zeros(len(self))

        for i,polarized_photon in enumerate(self):
            stokes[0, i] = polarized_photon.stokesVector().s0
            stokes[1, i] = polarized_photon.stokesVector().s1
            stokes[2, i] = polarized_photon.stokesVector().s2
            stokes[3, i] = polarized_photon.stokesVector().s3
            polarization_degrees[i] = polarized_photon.circularPolarizationDegree()

        array_dict["s0"] = stokes[0, :]
        array_dict["s1"] = stokes[1, :]
        array_dict["s2"] = stokes[2, :]
        array_dict["s3"] = stokes[3, :]
        array_dict["polarization degree"] = polarization_degrees

        return array_dict


    def toString(self):
        """Returns a string containing the parameters characterizing each polarized photon in the bunch."""
        bunch_string = str()

        for i in range(self.getNumberOfPhotons()):
            photon = self.getPhotonIndex(i)
            string_to_attach = str(photon.energy()) + " " + \
                               photon.unitDirectionVector().toString() + "\n"
            bunch_string += string_to_attach
        return bunch_string


#
#
#
class PolarizedPhotonBunchDecorator(PhotonBunchDecorator):

    def toDictionary(self):
        """Created a dictionary containing information about the bunch.

        Returns
        -------
        dict
            Information in tags: "number of photons", "energies", "deviations", "vx", "vy", "vz", "s0", "s1", "s2", "s3" and "polarization degree".

        """

        """defines a dictionary containing information about the bunch."""
        array_dict = PhotonBunch.toDictionary(self)

        stokes = numpy.zeros([4, len(self)])
        polarization_degrees = numpy.zeros(len(self))

        for i,polarized_photon in enumerate(self):
            stokes[0, i] = polarized_photon.stokesVector().s0
            stokes[1, i] = polarized_photon.stokesVector().s1
            stokes[2, i] = polarized_photon.stokesVector().s2
            stokes[3, i] = polarized_photon.stokesVector().s3
            polarization_degrees[i] = polarized_photon.circularPolarizationDegree()

        array_dict["s0"] = stokes[0, :]
        array_dict["s1"] = stokes[1, :]
        array_dict["s2"] = stokes[2, :]
        array_dict["s3"] = stokes[3, :]
        array_dict["polarization degree"] = polarization_degrees

        return array_dict

    def addPhoton(self, to_be_added):
        """Adds a photon to the bunch.

        Parameters
        ----------
        to_be_added : Photon instance

        """
        self.setEnergy(numpy.append(self.energy(), to_be_added.energy()))
        self.setUnitDirectionVector(self.unitDirectionVector().concatenate(to_be_added.unitDirectionVector()))
        self.setStokesVector(self.stokesVector().concatenate(to_be_added.stokesVector()))

    def getPhotonIndex(self, index):
        """Returns the photon in the bunch with a given index.

        Parameters
        ----------
        index : int
            The photon index to be referenced.

        Returns
        -------
        Photon instance
            The photon (referenced, not copied).

        """
        v = self.unitDirectionVector()
        vx = v.components()[0]
        vy = v.components()[1]
        vz = v.components()[2]
        s = self.stokesVector()
        return PolarizedPhoton(energy_in_ev=self.energy()[index],
                                      direction_vector=Vector(vx[index], vy[index], vz[index]),
                                      stokes_vector=StokesVector([s.s0[index], s.s1[index], s.s2[index], s.s3[index]])
                               )

    def setPhotonIndex(self, index, polarized_photon):
        """Sets the photon in the bunch with a given index.

        Parameters
        ----------
        index : int
            The photon index to be modified.

        polarized_photon : Photon instance
            The photon to be stored.

        """
        energy = self.energy()
        v = self.unitDirectionVector()
        vx = v.components()[0]
        vy = v.components()[1]
        vz = v.components()[2]
        s0 = self.stokesVector().s0
        s1 = self.stokesVector().s1
        s2 = self.stokesVector().s2
        s3 = self.stokesVector().s3

        energy[index] = polarized_photon.energy()
        vx[index] = polarized_photon.unitDirectionVector().components()[0]
        vy[index] = polarized_photon.unitDirectionVector().components()[1]
        vz[index] = polarized_photon.unitDirectionVector().components()[2]

        s0[index] = polarized_photon.stokesVector().s0
        s1[index] = polarized_photon.stokesVector().s1
        s2[index] = polarized_photon.stokesVector().s2
        s3[index] = polarized_photon.stokesVector().s3


        self.setEnergy(energy)
        self.setUnitDirectionVector(Vector(vx, vy, vz))
        self.setStokesVector(StokesVector([s0,s1,s2,s3]))


class PolarizedPhotonBunch(PolarizedPhoton, PolarizedPhotonBunchDecorator):
    """
    The PolarizedPhotonBunch is a PolarizedPhoton stack.

    It inheritates from PolarizedPhoton and uses stacks for more efficient stockage. Additional methods
    useful for stacks or bunches are defined in PolarizedPhotonBunchDecorator.

    Constructor.

    Parameters
    ----------
    polarized_photons : list, optional
        List of PolarizedPhoton instances.

    """
    def __init__(self, polarized_photons=None):

        if polarized_photons == None:
            super().__init__(energy_in_ev=[],
                             direction_vector=Vector([],[],[]),
                             stokes_vector=StokesVector([[],[],[],[]]))
        else:
            n = len(polarized_photons)
            if n == 0:
                super().__init__(energy_in_ev=[],
                                 direction_vector=Vector([], [], []),
                                 stokes_vector=StokesVector([[], [], [], []]))
            else:
                energy = numpy.zeros(n)
                s0 = numpy.zeros(n)
                s1 = numpy.zeros(n)
                s2 = numpy.zeros(n)
                s3 = numpy.zeros(n)
                for i,el in enumerate(polarized_photons):
                    energy[i] = el.energy()
                    vv = el.unitDirectionVector()
                    ss = el.stokesVector()
                    if i == 0:
                        v = Vector(
                            vv.components()[0],
                            vv.components()[1],
                            vv.components()[2],
                        )
                        s = StokesVector([ss.s0, ss.s1, ss.s2, ss.s3])
                    else:
                        v.append(vv)
                        s.append(ss)
                self.setEnergy(energy)
                self.setUnitDirectionVector(v)
                super().__init__(energy_in_ev=energy, direction_vector=v, stokes_vector=s)

    @classmethod
    def initializeFromPolarizedPhoton(cls, photon_stack):
        """Construct a polarized photon bunch from a polarized photon stack.

        Parameters
        ----------
        photon_stack : instance of PolarizedPhoton

        Returns
        -------
        PolarizedPhotonBunch instance

        """
        out = PolarizedPhotonBunch()
        out.setEnergy(photon_stack.energy())
        out.setUnitDirectionVector(photon_stack.unitDirectionVector())
        out.setStokesVector(photon_stack.stokesVector())
        return out

    @classmethod
    def initializeFromArrays(cls, energy=[], vx=[], vy=[], vz=[], s0=[], s1=[], s2=[], s3=[]):
        """Construct a polarized photon bunch from arrays with photon energies, directions and stokes components.

        Parameters
        ----------

        energies : list, numpy array
            the array with photon energy in eV.
        vx : list, numpy array
            the array with X component of the direction vector.
        vy : list, numpy array
            the array with Y component of the direction vector.
        vz : list, numpy array
            the array with Z component of the direction vector.
        s0 : list, numpy array
            the array with S0 stokes vector.
        s1 : list, numpy array
            the array with S1 stokes vector.
        s2 : list, numpy array
            the array with S2 stokes vector.
        s3 : list, numpy array
            the array with S4 stokes vector.

        Returns
        -------
        PolarizedPhotonBunch instance


        """
        bunch = PolarizedPhotonBunch()
        bunch.setEnergy(numpy.array(energy))
        bunch.setUnitDirectionVector(Vector(numpy.array(vx),
                                            numpy.array(vy),
                                            numpy.array(vz)))
        bunch.setStokesVector( StokesVector(element_list=
            [numpy.array(s0, dtype=float),
             numpy.array(s1, dtype=float),
             numpy.array(s2, dtype=float),
             numpy.array(s3, dtype=float)]))

        return bunch

if __name__ == "__main__":

    npoint = 10
    vx = numpy.zeros(npoint) + 0.0
    vy = numpy.zeros(npoint) + 1.0
    vz = numpy.zeros(npoint) + 0.1

    energy = numpy.zeros(npoint) + 3000.0

    photon_stack = PolarizedPhoton(energy, Vector(vx, vy, vz), StokesVector([vx,vx,vx,vx]))



    #
    # vector
    #


    photon_bunch3 = PolarizedPhotonBunch().initializeFromPolarizedPhoton(photon_stack)

    photon_bunch4 = PolarizedPhotonBunch().initializeFromArrays(
        energy=energy, vx=vx, vy=vy, vz=vz, s0=vx, s1=vx, s2=vx, s3=vx)

    #
    # check
    #

    print(">>>>>>>>>>>>>>>>>> 3")
    print(photon_bunch3.toDictionary()['s0'].shape)
    print(photon_bunch3.toDictionary())
    print(">>>>>>>>>>>>>>>>>> 4")
    print(photon_bunch4.toDictionary()['s0'].shape)
    print(photon_bunch4.toDictionary())