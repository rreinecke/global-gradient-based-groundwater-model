---
title: 'G³M-f a global gradient-based groundwater modelling framwork'
tags:
- groundwater
- gradient-based
- global
- conjugate-gradient
authors:
- name: Robert Reinecke
  orcid: 0000-0001-5699-8584
  affiliation: 1
affiliations:
- name: Goethe University Frankfurt
  index: 1

date: 02 January 2018
bibliography: paper.bib
---

# Summary
In order to represent groundwater-surface water interactions as well as the impact of capillary rise on evapotranspiration in global-scale hydrological models, it is necessary to simulate the location and temporal variation of the groundwater table.
This requires replacing simulation of groundwater dynamics using groundwater storage variations in individual grid cells (independent from the storage variation in neighbouring cells) with hydraulic head gradient-based groundwater modelling.

The global gradient-based groundwater model framework G³M-f is an extesible model framework.
Its main purpose is to be used as a main building block for the global groundwater mode G³M.
G³M is a newly developed gradient-based groundwater model which adapts MODFLOW [@harbaugh2005modflow] principles for the global scale.
It is written in C++ and intended to be coupled to the global hydrology model WaterGAP (http://watergap.de) [@alcamo2003development;@doll2003global;@doll2012impact;@doell2014global;@muller2014sensitivity], but can also be used for regional groundwater models and coupling to other hydrology models.
While it is intended to be used as a in memory coupled model it is also capable of running a standard standalone groundwater model.
The code is available on globalgroundwatermodel.org [@g3m].

# References
