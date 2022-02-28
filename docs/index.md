---
layout: default
title: Home
nav_order: 1
description: "Documentation of the global groundwater modeling framework"
permalink: /
---


[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

# Publications

[Global model description in GMD](https://www.geosci-model-dev.net/12/2401/2019/)

[Sensitivity Analysis in HESS](https://www.hydrol-earth-syst-sci.net/23/4561/2019/hess-23-4561-2019.html)


# Data and code dois
[![DOI](https://zenodo.org/badge/109667597.svg)](https://zenodo.org/badge/latestdoi/109667597)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1315471.svg)](https://doi.org/10.5281/zenodo.1315471)

# The global gradient-based groundwater model framework G³M-f
The global gradient-based groundwater model framework G³M-f is an extensible model framework.
Its main purpose is to be used as a main building block for the global groundwater mode G³M.
G³M is a newly developed gradient-based groundwater model which adapts MODFLOW [@harbaugh2005modflow] principles for the globalscale.
It is written in C++ and intended to be coupled to the global hydrology model WaterGAP (http://watergap.de), but can also be used for regional groundwater models and coupling to other hydrology models.
While it is intended to be used as a in memory coupled model it is also capable of running a standard standalone groundwater model.

TODO why this model is aswesome and which user groups.....

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes.

### Prerequisites

To compile the program, you will need:
```
clang >= 13 with openMP (currently gcc is not supported)
cmake >= 3.15.3
libboost >= 1.71
libGMP
libGtest
lcov
```
### Build
```
mkdir build
cd build
cmake ../
make
```

### Equations
The three dimensional flow of water through the porous material between the cells is solved as a partial differential equation.
Where K is the hydraulic conductivity [L/T] along the three axis, S the specific storage and W is the volumentric flux per unit volume in and out of the groundwater system.
The hydraulic conductivity between two cells is caluclated b yusing the harmonic mean.
The equation is solved using a conjugate gradient approach and an Incomplete LUT preconditioner.
![](https://latex.codecogs.com/gif.latex?\frac{\partial}{\partial&space;x}\left&space;(&space;K_{x}&space;\frac{\partial&space;h}{\partial&space;x}&space;\right&space;)&space;&plus;&space;\frac{\partial}{\partial&space;y}\left&space;(&space;K_{y}&space;\frac{\partial&space;h}{\partial&space;y}&space;\right&space;)&space;&plus;&space;\frac{\partial}{\partial&space;z}\left&space;(&space;K_{z}&space;\frac{\partial&space;h}{\partial&space;z}&space;\right&space;)&space;&plus;&space;W&space;=&space;S_{s}&space;\frac{\partial&space;h}{\partial&space;t} "Main equation")

Additonal information on the equations can be found in the very detailed MODFLOW documentation: [Modflow 2005](https://water.usgs.gov/ogw/modflow/MODFLOW.html)

### Boundary Conditions
G³M support multiple boundary condition types:
* No-flow boundary
* Static head boundary
* General head boundary
* Groundwater recharge
* Lakes
* Wetlands
* Different river approaches

New flows can be defined in Model/ExternalFlows.hpp.
The domain boundary is currently defined implicitly through the input grid as no-flow for grid files and as ocean boundary for irregual grids.
This behaviour can be changed in DataProcessing/Neighbouring.hpp.

## Built With

* [Eigen3](http://eigen.tuxfamily.org) - Doing the math magic
* [GTest](https://github.com/google/googletest) - Test framework
* [libboost](http://www.boost.org) - C++ magic
* [OpenMP](http://www.openmp.org) - Accelerator und Multi-Core support
* [GMP](https://gmplib.org) - Large numbers

## Contributing

Please read [CONTRIBUTING.md](CONTRIBUTING.md) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/your/project/tags). 

## Authors and Contributors


## Authors and Contributors

* **Robert Reinecke** <span id="badgeCont935"><script type="text/javascript" src="https://publons.com/mashlets?el=badgeCont935&rid=K-3693-2019&size=small"></script></span> - *Initial work* *Maintainer*
* **Daniel Kretschmer** - *Maintainer*
* **Sebastian Ackermann** - *Maintainer*

### Past Contributors

* **Alexander Wachholz** - *Documentation review*
* **Christoph Niemann** - *Spatial IDs* *Developer*

## License

This project is licensed under the GNU General Public License - see the [LICENSE](LICENSE) file for details.
Please note that the code contains a modified version of the Eigen3 library which is published under the [MPL 2.0](https://www.mozilla.org/en-US/MPL/2.0/).

## Acknowledgments

* [Modflow 2005](https://water.usgs.gov/ogw/modflow/MODFLOW.html) for their great documentation
* [Eigen3](http://eigen.tuxfamily.org) for their awesome framework
