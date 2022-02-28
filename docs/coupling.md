---
layout: default
title: Model Coupling
nav_order: 4
description: "How to use the model"
permalink: /coupling
---

# In memory coupling
G³M-f is written with the coupling to other models in mind.
In contrast to other model coupling efforts, it is not necessary to write out files in one model and read them in in another model.
You can diretly link G³M-f with your existing executable and by providing a class in your already existing model code that implements the gw_interface, you are free to call the simulate() function at any timestep you like.
Furthermore, the interface provides pointer containers and callbacks to transfer data in memory without the need to waste time on I/O.

Please contact us if you need advice.


