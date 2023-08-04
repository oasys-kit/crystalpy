
__authors__ = ["E Cappelli, M Glass, M Sanchez del Rio - ESRF ISDD Advanced Analysis and Modelling"]
__license__ = "MIT"
__date__ = "23/11/2016"

from setuptools import setup
#
# memorandum (for pypi)
#
# python setup.py sdist
# python -m twine upload dist/crystalpy....

#
# memorandum (for documentation with numpydoc style)
#
# install payment (https://stackoverflow.com/questions/24555327/how-can-i-produce-a-numpy-like-documentation)
# pyment -o numpydoc Vector.py  # apply to all files,
#   complete/edit by hand, see https://numpydoc.readthedocs.io/en/latest/format.html ...
# patch Vector.py < Vector.py.patch
# cd doc
# sphinx-quickstart  # needed only once...
# cd ..
# iterate:
# rm doc/crystalpy*.rst doc/modules.rst
# sphinx-apidoc -o doc crystalpy
# cd doc
# make clean html


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
