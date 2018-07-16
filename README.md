[![status](http://joss.theoj.org/papers/5fda5a279db561b6d4c597bbbe574867/status.svg)](http://joss.theoj.org/papers/5fda5a279db561b6d4c597bbbe574867)
[![DOI](https://zenodo.org/badge/109667597.svg)](https://zenodo.org/badge/latestdoi/109667597)
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![DOI:110.5194/gmd-2018-120](https://zenodo.org/badge/DOI/10.1007/978-3-319-76207-4_15.svg)](https://doi.org/10.5194/gmd-2018-120)


# The global gradient-based groundwater model framework G³M
The global gradient-based groundwater model framework G³M-f is an extensible model framework.
Its main purpose is to be used as a main building block for the global groundwater mode G³M.
G³M is a newly developed gradient-based groundwater model which adapts MODFLOW [@harbaugh2005modflow] principles for the globalscale.
It is written in C++ and intended to be coupled to the global hydrology model WaterGAP (http://watergap.de), but can also be used for regional groundwater models and coupling to other hydrology models.
While it is intended to be used as a in memory coupled model it is also capable of running a standard standalone groundwater model.

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes.

### Prerequisites

```
clang >= 3.8 with openMP (currently gcc is not supported)
libboost >= 1.56
libGMP
libGtest
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

## Quick start
The following picture shows the conceptual example model:
![](docs/simple_model.png)

After compilation run:
```
simple_model
```
It will yield a depth to water table CSV file called wtd.csv for a simple model.

## How to use
The follwing will guide you through the building blocks of the simple model shipped along with the code.
It assumes that you've constructed your model domain and have input data for the following
* Groundwater recharge (recharge_simple.csv)
* Surface elevation (elevation_simple.csv)
* Rivers, location, elevation and depth (rivers_simple.csv)
* Hydrogeology (lithology_simple.csv)
* Riverbed conductance (rivers_simple.csv)
* Inital head guess (otherwise the model assumes the surface elevation as best guess) ("heads_simple.csv")

Center building stone for the framework is the GW_interface connecting any model with the groundwater code.
Implement this interface if you want to couple your model to G³M-f or build a custom standalone application.
In tests/SimpleModel you'll find an example implementation explained further in the following.

```
class GW_Interface {
    public:
        virtual ~GW_Interface() {}

        virtual void
        loadSettings() = 0;

        virtual void
        setupSimulation() = 0;

        virtual void
        writeData() = 0;

        virtual void
        simulate() = 0;
};
```

The following shows the code for a simple model loop running a steady-state model with daily timesteps.
```
void StandaloneRunner::simulate() {
    Simulation::Stepper stepper = Simulation::Stepper(_eq, Simulation::DAY, 1);
    for (Simulation::step step : stepper) {
        LOG(userinfo) << "Running a steady state step";
        step.first->toogleSteadyState();
        step.first->solve();
        sim.printMassBalances();
    }
    DataProcessing::DataOutput::OutputManager("data/out_simple.json", sim).write();
    //sim.save();
}
```

### Write out data
Write out of data is specified by a JSON file called out.json.
If you want to add custom fields you can do so in src/DataProcessing/DataOutput.
```
{
  "output": {
    "StaticResult": [
      {
        "name": "wtd",
        "type": "csv",
        "field": "DepthToWaterTable",
        "ID": "false",
        "position": "true"
      }
    ],
    "InnerIteration": {
    },
    "OuterIteration": {
    }
  }
}

```

### Config model
In order to configure the model variables you can simply change the .json file. Allowing you to change the convergence criteria and the location for your input files.

### Parameters
The following explains the main config parameters.

* model_config
  * nodes: A file describing the input grid
  * row_cols: true: neighbouring is determined by their position in an evenly grid, false: neighbouring is determined by their lat and lon position (currently only supports 5' resolution)
  * threads: How many computation threads to use to solve the equation
  * layers: Number of layers of the model domain
  * confinement: Define which of the layers is a confined layer
* numerics
  * solver: Currently only Preconditioned Conjugent Gradient, code for a newton approach is available but untested
  * iterations: Number of picard iterations
  * closingcrit: Inf norm of the residuals
  * headchange: Closing criterion for max. head change for 3 consecutive iterations
  * damping: Damping of residuals in between picard iterations
* input: Internaly the model code assumes time dependant parameters to be per day
  * data_config: Describes wether default data is used or a input file should be read
  * default_data: specifiy default parameters
  * data: Inputdata - can be modified according to the users need. The shown inputs are the supported defaults

```
{
  "config": {
    "model_config": {
      "nodes": "grid_simple.csv",
      "row_cols": "true",
      "stadystate": "true",
      "numberofnodes": 100,
      "threads": 1,
      "layers": 2,
      "confinement": [
        "false",
        "true"
      ],
      "cache": "false",
      "adaptivestepsize": "false",
      "boundarycondition": "SeaLevel",
      "sensitivity": "false"
    },
    "numerics": {
      "solver": "PCG",
      "iterations": 500,
      "inner_itter": 10,
      "closingcrit": 1e-8,
      "headchange": 0.0001,
      "damping": "false",
      "min_damp": 0.01,
      "max_damp": 0.5,
      "stepsize": "daily"
    },
  "input": {
    "data_config": {
      "k_from_lith": "true",
      "k_ocean_from_file": "false",
      "specificstorage_from_file": "false",
      "specificyield_from_file": "false",
      "k_river_from_file": "true",
      "aquifer_depth_from_file": "false",
      "initial_head_from_file": "true",
      "data_as_array": "false"
    },
    "default_data": {
      "initial_head": 5,
      "K": 0.008,
      "oceanK": 800,
      "aquifer_thickness": [
        10,
        10
      ],
      "anisotropy": 10,
      "specificyield": 0.15,
      "specificstorage": 0.000015
    },
    "data": {
      "recharge": "recharge_simple.csv",
      "elevation": "elevation_simple.csv",
      "rivers": "rivers_simple.csv",
      "lithologie": "lithology_simple.csv",
      "river_conductance": "rivers_simple.csv",
      "initial_head": "heads_simple.csv"
    }
  }
  }
}
```

## Deployment in other models
The main steps towards your own model is to implement the GW_interface and provide a DataReader.
A standlone version can be easily implemented by extending the simple example provided above.

### In memory coupling
G³M-f is written with the coupling to other models in mind.
In contrast to other model coupling efforts, it is not necessary to write out files in one model and read them in in another model.
You can diretly link G³M-f with your existing executable and by providing a class in your already existing model code that implements the gw_interface, you are free to call the simulate() function at any timestep you like.
Furthermore, the interface provides pointer containers and callbacks to transfer data in memory without the need to waste time on I/O.

Please contact us if you need advice.

## Running the tests
Automated tests consits of gunit test which are compiled automatically with the attached cmake file.
You can run them by executing the test executable.

```
runUnitTests
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

* **Robert Reinecke** - *Initial work*

## License

This project is licensed under the GNU General Public License - see the [LICENSE](LICENSE) file for details.
Please note that the code contains a modified version of the Eigen3 library which is published under the [MPL 2.0](https://www.mozilla.org/en-US/MPL/2.0/).

## Acknowledgments

* [Modflow 2005](https://water.usgs.gov/ogw/modflow/MODFLOW.html) for their great documentation
* [Eigen3](http://eigen.tuxfamily.org) for their awesome framework
