---
layout: default
title: Data input
nav_order: 2
parent: Input and Outputs
---

# Data Input

## Model grid (e.g., grid.csv)
Model grid definition. X/Y give the lon/lat position. The remaining data sets will read only data for spatIDs read from this file.  

[] 1st column: "spatID" - spatial ID

[] 2nd column: "X" - longitude

[] 3rd column: "Y" - latitude

[] 4th column: "area" - area

[] (5th column: "col" - ID of column in grid - optional for small models)

[] (6th column: "row" - ID of row in grid - optional for small models)

```
spatID,X,Y,area,col,row
2247710,172.957001,-34.375099,70.552101,78,0
2247711,173.039993,-34.375099,70.552101,79,0
2248179,172.957001,-34.458401,70.481796,78,1
...
```

## Groundwater recharge (e.g., recharge.csv)
Groundwater recharge rate.

[] 1st column: "spatID" - spatial ID

[] 2nd column: "data" - groundwater recharge in meters per day

```
spatID,data
2247710,0.018
2247711,0.02
2248179,0.002
...
```

## Surface elevation (e.g., elevation.csv)
Elevation of the earths surface.

[] 1st column: "spatID" - spatial ID

[] 2nd column: "data" - surface elevation in meters

```
spatID,data
2247710,0.018
2247711,0.02
2248179,0.002
...
```

## Surface water elevation (e.g., elevation_30.csv)
Elevation of surface waters: rivers, wetlands, lakes. The 30th percentile of the elevation range with in respective grid cell on a finer resolution gave reasonable results for the global model. 

[] 1st column: "spatID" - spatial ID

[] 2nd column: "data" - surface water elevation in meters

```
spatID,data
2247710,0.018
2247711,0.02
2248179,0.002
...
```

## Rivers (e.g., rivers.csv)

[] 1st column: "spatID" - spatial ID

[] 2nd column: "Head" - river head in meters

[] 3rd column: "Bottom" - river bottom elevation in meters

[] 4th column: "Conduct" - riverbed conductance in meters per day

```
spatID,Head,Bottom,Conduct
2247710,3,-10,48
2247711,3,-10,48
2248179,3,-10,48
...
```

## Permeability (e.g., lithology.csv)
Permeability of the aquifer.

[] 1st column: "spatID" - spatial ID

[] 2nd column: "data" - permeability in meters per day

```
spatID,data
2247710,0.5
2247711,2
2248179,1
...
```

## Initial heads (otherwise the model assumes the surface elevation as best guess) (e.g., initial_heads.csv)
Initial head of the groundwater. May be a first guess and/or based on observations. 

[] 1st column: "spatID" - spatial ID

[] 2nd column: "data" - initial head in meters

```
spatID,data
2247710,-10
2247711,-12
2248179,-30
...
```

## Water table depth (e.g., water_table_depth.csv)
Initial water table depth of the groundwater. May be a first guess and/or based on observations. Either water table depth OR initial heads should be used as input data. Water table depth and elevation are used to compute the initial heads.    

[] 1st column: "spatID" - spatial ID

[] 2nd column: "data" - water table depth in meters

```
spatID,data
2247710,-10
2247711,-12
2248179,-30
...
```

## Slope (e.g., slope.csv)
Terrain slope.

[] 1st column: "spatID" - spatial ID

[] 2nd column: "data" - slope in degrees

```
spatID,data
2247710,0.104
2247711,0.186
2248179,0.153
...
```

## E-Folding (e.g., efolding.csv)
E-folding factor f used by Fan et al. (2013)[https://www.science.org/doi/10.1126/science.1229881] to calculate the conductivity of lower layers by multiplying the upper layer value by exp(-50m f^-1)^-1.

[] 1st column: "spatID" - spatial ID

[] 2nd column: "data" - e-folding factor 

```
spatID,data
2247710,62.413
2247711,60.538
2248179,41.531
...
```

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
