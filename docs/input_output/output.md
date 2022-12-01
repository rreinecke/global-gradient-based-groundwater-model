---
layout: default
title: Outputs
nav_order: 10
parent: Input and Outputs 
---

# Write out data
The settings in tests/out.json define the export of results:
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
