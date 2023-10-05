__authors__ = ["M Sanchez del Rio, E Cappelli, M Glass  - ESRF ISDD MEG Modelling group"]
__license__ = "MIT"
__date__ = "2016-2023"

from setuptools import setup

#
# memorandum (for developer installation)
#
# git clone https://github.com/oasys-kit/crystalpy
# cd crystalpy
# python -m pip install -e . --no-deps --no-binary :all:

#
# memorandum (for updating pypi)
#
# python setup.py sdist
# python -m twine upload dist/crystalpy....

#
# memorandum (for documentation with numpydoc style)
#
# ** install sphinx:
# pip install sphinx
# pip install sphinxcontrib-apidoc
# pip install sphinx-rtd-theme
# pip install nbsphinx
# install payment (https://stackoverflow.com/questions/24555327/how-can-i-produce-a-numpy-like-documentation)
# ** create embeded doc
# pyment -o numpydoc Vector.py  # apply to all files,
#   complete/edit by hand, see https://numpydoc.readthedocs.io/en/latest/format.html ...
# patch Vector.py < Vector.py.patch

# ** some initialization
# sphinx-quickstart  # needed only once...
# ** iterate
# rm docs/crystalpy*.rst docs/modules.rst
# sphinx-apidoc -o docs crystalpy
# make clean html


setup(name='crystalpy',
      version='0.0.17',
      description='Python crystal polarization calculation',
      author='Manuel Sanchez del Rio, Edoardo Cappelli, Mark Glass',
      author_email='srio@esrf.eu',
      url='https://github.com/oasys-kit/crystalpy/',
      packages=['crystalpy',
                'crystalpy.util',
                'crystalpy.diffraction',
                'crystalpy.polarization'],
      install_requires=[
                        'numpy',
                        'scipy',
                        'mpmath',
                        'dabax',
                       ],
      test_suite='tests'
      )
