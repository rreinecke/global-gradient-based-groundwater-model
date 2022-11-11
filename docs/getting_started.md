---
layout: default
title: Getting Started
nav_order: 3
has_children: true
description: "How to use the model"
permalink: /getstarted
---

# Prerequisites

To compile the program, you will need:
```
clang >= 13 with openMP (currently gcc is not supported)
cmake >= 3.15.3
libboost >= 1.71
libGMP
libGtest
lcov
```

# How to compile
```
mkdir build
cd build
cmake ../
make
```


Check out [here](http://globalgroundwatermodel.org/tutorial/tut1) how to set up a simple model.




