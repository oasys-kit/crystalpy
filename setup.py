
__authors__ = ["E Cappelli, M Glass, M Sanchez del Rio - ESRF ISDD Advanced Analysis and Modelling"]
__license__ = "MIT"
__date__ = "23/11/2016"

from setuptools import setup
#
# memorandum (for pypi)
#
# python setup.py sdist upload



setup(name='crystalpy',
      version='0.0.11',
      description='Python crystal polarization calcution',
      author='Manuel Sanchez del Rio, Edoardo Cappelli, Mark Glass',
      author_email='srio@esrf.eu',
      url='https://github.com/oasys-kit/crystalpy/',
      packages=['crystalpy',
                'crystalpy.util',
                'crystalpy.diffraction',
                'crystalpy.polarization',
                'crystalpy.examples',
                'crystalpy.tests'],
      install_requires=[
                        'numpy',
                        'scipy',
                        'mpmath',
                        'dabax'
                       ],
      test_suite='tests'
      )
