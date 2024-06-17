---
layout: default
title: Coupling with GHMs
nav_order: 4
description: "How to use the model"
permalink: /coupling
---
## Deployment in other models
The main steps towards coupling G³M-f with your own model is to (1) implement the GW_interface and (2) provide a DataReader.

# In memory coupling
G³M-f is written with the coupling to other models in mind.
In contrast to other model coupling efforts, it is not necessary to write out files in one model and read them in in another model.
You can diretly link G³M-f with your existing executable and by providing a class in your already existing model code that implements the gw_interface, you are free to call the simulate() function at any timestep you like.
Furthermore, the interface provides pointer containers and callbacks to transfer data in memory without the need to waste time on I/O.

[Here](http://globalgroundwatermodel.org/getting_started/coupling) you can find a tutorial on how to couple. 

Please contact us if you need advice.

