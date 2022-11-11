---
layout: default
title: Home
nav_order: 1
description: "Documentation of the global groundwater modeling framework."
permalink: /
---
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

# The framework G³M-f
The global gradient-based groundwater model framework (G³M-f) is an extensible program to build groundwater models.
Its primary purpose is to be used as a central building block for the global groundwater model G³M.
G³M is a global gradient-based groundwater model which adapts MODFLOW principles for the global scale.
It is written in C++ and intended to be coupled to the global hydrology model WaterGAP (http://watergap.de).
While it is intended to be used as a memory coupled model, it can also run a standard standalone groundwater model.

It is mainly intended to be used by developers of global hydrological models (including land surface models and others) who want to simulate groundwater in a gradient-based manner.
In principle, it can also implement any other, e.g., regional, groundwater model.
Its main feature compared to MODFLOW is its extensibility, speed, and in-memory coupling capability.

For collaboration on coupling this model to other hydrological models, please contact Dr. Robert Reinecke (rreinecke on github).

# Publications

[Global model description in GMD](https://www.geosci-model-dev.net/12/2401/2019/)

[Sensitivity Analysis in HESS](https://www.hydrol-earth-syst-sci.net/23/4561/2019/hess-23-4561-2019.html)


# Data and code dois
Code publication in JOSS: [![DOI](https://zenodo.org/badge/109667597.svg)](https://zenodo.org/badge/latestdoi/109667597)

Outputs of Reinecke et al 2019 in GMD: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1315471.svg)](https://doi.org/10.5281/zenodo.1315471)

# Check out the project on GitHub
[GitHub](https://github.com/rreinecke/global-gradient-based-groundwater-model)


# Prerequisites

To compile the program, you will need:
```
clang >= 13 with openMP (currently gcc is not supported)
cmake >= 3.15.3
libboost >= 1.71
libGMP
libGtest
lcov
```



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
* **Daniel Kretschmer** - *Variable density flow* *Maintainer*
* **Sebastian Ackermann** - *WaterGAP coupling* *Maintainer*

### Past Contributors

* **Alexander Wachholz** - *Documentation review and NZ model* *Developer*
* **Christoph Niemann** - *Spatial IDs* *Developer*

## License

This project is licensed under the GNU General Public License - see the [LICENSE](LICENSE) file for details.
Please note that the code contains a modified version of the Eigen3 library which is published under the [MPL 2.0](https://www.mozilla.org/en-US/MPL/2.0/).

## Acknowledgments

* [Modflow 2005](https://water.usgs.gov/ogw/modflow/MODFLOW.html) for their great documentation
* [Eigen3](http://eigen.tuxfamily.org) for their awesome framework
