---
layout: default
title: Model Concept
nav_order: 2
has_children: true
description: ""
permalink: /concept
---

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

