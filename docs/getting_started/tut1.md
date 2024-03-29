---
layout: default
title: "Tutorial 1 - a simple model"
parent: Getting Started
nav_order: 1
---

# Tutorial 1 - a simple model

## simple conceptual model

The following picture shows the conceptual example model:
![](simple_model.png)

## Running the model

After compilation (see [here](http://globalgroundwatermodel.org/getstarted)):
```
simple_model
```
It will yield a depth to water table CSV file called wtd.csv for a simple model. The output can be changed (see "Write out data" below)

## How to use
The following will guide you through the building blocks of the simple model shipped along with the code.
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
Data is written into a file depending on the settings in build/data/out_simple.json:
* "name": name of the created  output file
* "type": file type (current options are "csv" and "gfs-json", help for netCDF output implementation is appreciated)
* "field": name of the field (the list of field options can be found in src/DataProcessing/DataOutput/FieldCollector.hpp, e.g., "Velocity")
* "ID": if true, the node ID is written into the output file
* "position": if true, node positions (Y and X, or latitude and longitude) are written into the file

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
The output file (here "wtd.csv") is written into the build directory.

For advanced users: If you want to add custom fields you can do so in src/DataProcessing/DataOutput.

### Model configuration

The model parameters (e.g. aquifer settings, convergence criteria, location of input files) can be configured by changing tests/SimpleModel/config.json.

#### Main config parameters

* model_config
  * nodes: A file describing the input grid
  * row_cols: if true, neighbouring is determined by their position in an even grid; if false, neighbouring is determined by their lat and lon position (currently only supports 5' resolution)
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

#### The config.json file
```
{
  "config": {
    "model_config": {
      "nodes": "grid_simple.csv",
      "row_cols": "true",
      "steadystate": "true",
      "number_of_nodes": 100,
      "number_of_rows": 10,
      "number_of_cols": 10,
      "edge_length_rows": 3.162277,
      "edge_length_cols":3.162277,
      "threads": 1,
      "layers": 2,
      "one_layer_approach": "false",
      "confinement": [
        "true",
        "true"
      ],
      "cache": "false",
      "adaptivestepsize": "false",
      "boundarycondition": "GeneralHeadBoundary",
      "sensitivity": "false"
    },
    "numerics": {
      "solver": "PCG",
      "iterations": 100,
      "inner_itter": 10,
      "closingcrit": 1e-90,
      "headchange": 0.000001,
      "damping": "false",
      "min_damp": 0.01,
      "max_damp": 0.5,
      "stepsize": "daily",
      "wetting_approach": "nwt"
    },
  "input": {
    "data_config": {
      "k_from_lith": "true",
      "k_ghb_from_file": "false",
      "specificstorage_from_file": "false",
      "specificyield_from_file": "false",
      "k_river_from_file": "true",
      "aquifer_depth_from_file": "false",
      "initial_head_from_file": "true",
      "data_as_array": "false"
    },
    "default_data": {
      "initial_head": 100,
      "K": 0.008,
      "ghb_K": 800,
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
      "lithology": "lithology_simple.csv",
      "river_conductance": "rivers_simple.csv",
      "initial_head": "heads_simple.csv"
    }
  }
  }
}
```

## Further documentation:

[The full documentation](html_doku/index.html)
