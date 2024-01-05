.. currentmodule:: crystalpy

===========
Package API
===========
This page lists main classes in this package.


``crystalpy`` classes.


diffraction
-----------
``crystalpy.diffraction`` classes for diffraction setup and calculations

.. autosummary::
   :toctree: generated/

   diffraction.Diffraction
   diffraction.DiffractionSetupAbstract
   diffraction.DiffractionSetupSweeps
   diffraction.DiffractionSetupDabax
   diffraction.DiffractionSetupXraylib
   diffraction.DiffractionSetupShadowPreprocessorV1
   diffraction.DiffractionSetupShadowPreprocessorV2
   diffraction.DiffractionExceptions
   diffraction.DiffractionResult
   diffraction.GeometryType
   diffraction.PerfectCrystalDiffraction

polarization
------------
``crystalpy.polarization`` classes for polarization calculations

.. autosummary::
   :toctree: generated/

   polarization.CrystalPhasePlate
   polarization.MuellerDiffraction
   polarization.MuellerMatrix
   polarization.MuellerResult

util
----
``crystalpy.util`` classes and methods with utilities and tools

.. autosummary::
   :toctree: generated/

   util.Vector
   util.StokesVector
   util.Photon
   util.PhotonBunch
   util.PolarizedPhoton
   util.PolarizedPhotonBunch
   util.ComplexAmplitudePhoton
   util.ComplexAmplitudePhotonBunch
   util.CalculationStrategy
   util.bragg_preprocessor_file_io
