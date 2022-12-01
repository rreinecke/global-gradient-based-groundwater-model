---
layout: default
title: Data input
nav_order: 2
parent: Input and Outputs
---

# Data Input

## Model grid (grid.csv)
1st column: "spatID" - spatial ID
2nd column: "X" - longitude
3rd column: "Y" - latitude
4th column: "area" - area
(5th column: "col" - grid column - optional for small models)
(6th column: "row" - grid row - optional for small models)

spatID,X,Y,area,col,row
2247710,172.957001,-34.375099,70.552101,78,0
2247711,173.039993,-34.375099,70.552101,79,0
2248179,172.957001,-34.458401,70.481796,78,1

## Groundwater recharge (recharge.csv)
1st column: "spatID" - spatial ID
2nd column: "data" - groundwater recharge in meters per day

spatID,data
2250266,0.574833277
2249872,0.574833277
2249873,0.574833277

## Surface elevation (elevation.csv)
1st column: spatial ID
2nd column: surface elevation in meters

## Rivers, location, elevation and depth (rivers.csv)


## Permeability (lithology.csv)

## Riverbed conductance (rivers.csv)

## Inital heads guess (otherwise the model assumes the surface elevation as best guess) ("initial_heads.csv")
