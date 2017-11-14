# The global gradient-based groundwater model framework G³M
The global gradient-based groundwater model framework G³M-f is an extesible model framework that is the basis for the G³M coupled to the global hydrologic model WaterGAP (http://watergap.de/).

While it is intended to be used as a in memory coupled model it is also capable of running a standard standalone groundwater model.

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes.

### Prerequisites

```
gcc > 6.2 or clang > 3.8 with openMP
libboost >= 1.56
OpenMP
libGMP
Gtest
```

### How to use
Center building stone for the framework is the GW_interface connecting any model with the groundwater code.

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

## Reading in Data

```
TODO
```

## Write out data
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

## Config model
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

## Running a simple model
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
## Deployment in other models
Just implement the GW_interface and provide a DataReader.

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

Please read [CONTRIBUTING.md]() for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/your/project/tags). 

## Authors

* **Robert Reinecke** - *Initial work*

See also the list of [contributors]() who participated in this project.

## License

This project is licensed under the GNU General Public License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* [Modflow 2005](https://water.usgs.gov/ogw/modflow/MODFLOW.html) for their great documentation
* [Eigen3](http://eigen.tuxfamily.org) for their awesome framework
