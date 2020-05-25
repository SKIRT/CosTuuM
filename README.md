[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3842422.svg)](https://doi.org/10.5281/zenodo.3842422)

# C++ T-Matrix code
Based on the T-matrix code by Mishchenko (https://www.giss.nasa.gov/staff/mmishchenko/t_matrix.html).

© Astronomical Observatory, Ghent University

The code and documentation in this repository is open source and freely available to the worldwide scientific community (see the [license file](COPYING) in this repository). We invite you to use CosTuuM or other portions of this project in your work. If you do so, we kindly ask that you cite the relevant papers.

We also explicitly welcome your contributions to the project, as set out in our [Guidelines for Contributing](CONTRIBUTING.md). A contribution can take many forms: asking a question, reporting a bug, suggesting a new feature, correcting or improving the documentation, writing a tutorial, fixing a bug, implementing a new class or a new module, and so on. Please adhere to our [Code of Conduct](CODE_OF_CONDUCT.md) to help maintain a constructive and open environment for all.

## Third party libraries

CosTuuM includes code from the following third-party libraries:
 - [QuickSched](https://gitlab.cosma.dur.ac.uk/swift/quicksched), a scheduler for task-based shared-memory parallelism ([Gonnet et al., 2016](https://arxiv.org/abs/1601.05384)), released under the GNU General Public License v3.0. The original license is included in [the relevant subdirectory](quicksched/COPYING). Small changes were made to QuickSched to make it compatible with C++; these changes have been indicated in the relevant files.
 - [CMacIonize](https://github.com/bwvdnbro/CMacIonize), a Monte Carlo photoionization and radiation hydrodynamics code ([Vandenbroucke & Wood, 2018](https://ui.adsabs.harvard.edu/abs/2018ascl.soft02003V/abstract)), released under the same license as CosTuuM. Various files were copied from CMacIonize without changes; this has been indicated in these files.

## Quick start guide

The code can be configured using CMake:

```
cmake <PATH TO SOURCE FOLDER>
```

we recommend an out-of-source configuration (where the build folder is 
either a subdirectory of the source folder or a folder that is outside 
the source folder). Once the configuration has succeeded, the code can 
be compiled using Make:

```
make -j <NUMBER OF THREADS TO USE>
```

The CosTuuM Python library can be installed for the current user using

```
make install
```

This will make it available from within Python, simply

```
import CosTuuM
```
