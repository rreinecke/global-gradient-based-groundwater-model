---
layout: default
title: Input - Configuration file
nav_order: 1
parent: Input and Outputs
---


# Model configuration

The model parameters (e.g. aquifer settings, convergence criteria, location of input files) can be configured by changing the config.json (an example is shown below).

## Main config parameters

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

## Example config.json file 
This is the same as in Tutorial 1 - a simple model.
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
