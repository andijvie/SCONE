# SCONE FORK FOR ANDY QITIAN ZHANG MPHIL DISSERTATION: Fission Matrix Acceleration in Monte Carlo Simulations and Its Impact on Neutron Clustering
This version is forked from the 6 May 2026 revision of the SCONE main branch, available at “CambridgeNuclear” on GitHub.

## List of changes (exhaustive):
 * FISSION MATRIX ACCELERATION: has been implemented for inactive and active cycles. The neutron weights are scaled to mach the fundamental EV of the FM. Tallies are implicit: based on the expected number of neutrons created from every collision.
 * FISSION MATRIX ACCELERATION: can use cumulative tally (default) or set moving window size
 * FISSION MATRIX ACCELERATION: setting to skip tallying for some initial number of cycles
 * FISSION MATRIX ACCELERATION: can force the eigenvector to always be homogeneous (for testing)
 * FISSION MATRIX ACCELERATION: has a setting for debug messages that show the calculations.
 * FISSION MATRIX ACCELERATION: can set a minimum number of cycles tallied (in eigenphysics package)
 * NEW TALLY CLERK: tallies fission matrix (need to enable as a setting)
 * NEW TALLY CLERK: collisionClerkCycle prints the (cycle or cumulative) collision flux in every MC cycle.
 * NEW EIGENPHYSICS SETTING: preserve total weight? (automatically enabled whith FM acceleration)
 * NEW EIGENPHYSICS SETTING: use combing as normalisation method? (automatically enabled whith FM acceleration)
 * NEW EIGENPHYSICS SETTING: print the FM EV of every cycle to a file?
 * NEW EIGENPHYSICS SETTING: print the fission source to a file?
 * NEW EIGENPHYSICS SETTING: show debug messages?
 * The output file can now handle magnitudes E+-999
 * Addition of 1D Fissile-Nonfissile-Fissile system stored as FNFslab (uses absorb_XSS and fissile_XSS)
 * Addition of 1D homogeneous system stored as FNFslab (uses homg_XSS)


# SCONE
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)
[![build-and-test-ubuntu](https://github.com/CambridgeNuclear/SCONE/actions/workflows/build-and-test.yml/badge.svg)](https://github.com/CambridgeNuclear/SCONE/actions/workflows/build-and-test.yml)
[![Documentation Status](https://readthedocs.org/projects/scone/badge/?version=latest)](https://scone.readthedocs.io/en/latest/?badge=latest)

SCONE (**S**tochastic **C**alculator **O**f **N**eutron Transport **E**quation) is an object-oriented Monte Carlo
particle transport code for reactor physics. It is intended as an accessible environment for
graduate students to test and develop their ideas before contributing them to more established
codes suitable for design calculations.

SCONE documentation is hosted at: <https://scone.readthedocs.io>

To cite SCONE, please use the following:
```bibtex
@article{sconeANE,
  title   = {Status of the SCONE Monte Carlo neutron transport code},
  journal = {Annals of Nuclear Energy},
  volume  = {227},
  pages   = {112015},
  year    = {2026},
  issn    = {0306-4549},
  doi     = {10.1016/j.anucene.2025.112015},
  url     = {https://www.sciencedirect.com/science/article/pii/S0306454925008321},
  author  = {Valeria Raffuzzi and Paul Cosgrove and Mikolaj Adam Kowalski},
  keywords = {SCONE, Monte Carlo, Neutron transport},
}
```

## Prerequisites
Required

* Cmake (>=3.10)
* Fortran compiler, gfortran (>=8.3)
* LAPACK and BLAS Libraries
* GNU/Linux operating system

Optional

* pFUnit 4 test framework
* Python 3 interpreter

## Installation
Instructions are available in the Sphinx documentation.

## Compiling Documentation
Sphinx documentation is available in the docs folder. It is readable with any reStructuredText (RST)
viewer, but it is best to compile to html.

Compiling documentation requires few python packages. You can install them all with the following
command. Option `--user` installs them in your home directory and does not require administrator access.
```
pip install --user -U sphinx, sphinx_rtd_theme
```
Then navigate to `docs` folder and compile using `make`
```
make html
```

HTML documentation should now be available in `./_build/html`

## Licence
This project is licensed under MIT Licence - see the [LICENCE](LICENCE) file for details.
