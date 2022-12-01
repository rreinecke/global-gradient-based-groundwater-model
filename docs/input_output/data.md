---
layout: default
title: Data input
nav_order: 2
parent: Input and Outputs
---

# Data Input

## Model grid (e.g., grid.csv)
1st column: "spatID" - spatial ID
2nd column: "X" - longitude
3rd column: "Y" - latitude
4th column: "area" - area
(5th column: "col" - grid column - optional for small models)
(6th column: "row" - grid row - optional for small models)

```
spatID,X,Y,area,col,row
2247710,172.957001,-34.375099,70.552101,78,0
2247711,173.039993,-34.375099,70.552101,79,0
2248179,172.957001,-34.458401,70.481796,78,1
...
```

## Groundwater recharge (e.g., recharge.csv)
1st column: "spatID" - spatial ID
2nd column: "data" - groundwater recharge in meters per day

```
spatID,data
2247710,0.018
2247711,0.02
2248179,0.002
...
```

## Surface elevation (e.g., elevation.csv)
1st column: "spatID" - spatial ID
2nd column: "data" - surface elevation in meters

```
spatID,data
2247710,0.018
2247711,0.02
2248179,0.002
...
```

## Surface water elevation (e.g., elevation_30.csv)
1st column: "spatID" - spatial ID
2nd column: "data" - surface water elevation in meters (e.g., the 30th percentile of the elevation range on a finer resolution)

```
spatID,data
2247710,0.018
2247711,0.02
2248179,0.002
...
```

## Rivers (e.g., rivers.csv)
1st column: "spatID" - spatial ID
2nd column: "Head" - river head in meters
3rd column: "Bottom" - river bottom elevation in meters
4th column: "Conduct" - riverbed conductance in meters per day

```
spatID,Head,Bottom,Conduct
2247710,3,-10,48
2247711,3,-10,48
2248179,3,-10,48
...
```

## Permeability (e.g., lithology.csv)
1st column: "spatID" - spatial ID
2nd column: "data" - permeability in meters per day

```
spatID,data
2247710,0.5
2247711,2
2248179,1
...
```

## Inital heads guess (otherwise the model assumes the surface elevation as best guess) (e.g., initial_heads.csv)
1st column: "spatID" - spatial ID
2nd column: "data" - initial head in meters

```
spatID,data
2247710,-10
2247711,-12
2248179,-30
...
```

## Water table depth (e.g., water_table_depth.csv)
1st column: "spatID" - spatial ID
2nd column: "data" - water table depth in meters

```
spatID,data
2247710,-10
2247711,-12
2248179,-30
...
```

## Slope (e.g., slope.csv)
1st column: "spatID" - spatial ID
2nd column: "data" - slope in degrees

```
spatID,data
2247710,0.104
2247711,0.186
2248179,0.153
...
```

## E-Folding (e.g., efolding.csv)





## CMakeLists for input data
To ensure the right data is used in the current model, the CMakeLists file inside the model folder needs to contain a line copying each desired file to the data folder of the build directory:
```
configure_file(input_file_name.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/input_file_name.csv COPYONLY)
```

The input data may be stored where it suits the current projects' purpose. Just adapt the CMakeLists file.

In the case of tutorial 1, the respective CMakeLists file is in the folder tests/SimpleModel/. These are the lines that define where to copy the data to:

```
configure_file(data/elevation.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/elevation_simple.csv COPYONLY)
configure_file(data/grid.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/grid_simple.csv COPYONLY)
configure_file(data/heads.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/heads_simple.csv COPYONLY)
configure_file(data/lithology.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/lithology_simple.csv COPYONLY)
configure_file(data/recharge.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/recharge_simple.csv COPYONLY)
configure_file(data/rivers.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/rivers_simple.csv COPYONLY)
```


## Further data
Further data can be used as an input, such as the slope, mean water table depth and information about aquifer folding.
