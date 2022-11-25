set(SOURCE
   ${SOURCE}
   ${CMAKE_CURRENT_SOURCE_DIR}/runner.cpp
   PARENT_SCOPE
)
set(HEADERS
   ${HEADERS}
   ${CMAKE_CURRENT_SOURCE_DIR}/runner.hpp
   ${CMAKE_CURRENT_SOURCE_DIR}/GlobalDataReader.hpp
   PARENT_SCOPE
)
configure_file(config.json ${CMAKE_CURRENT_BINARY_DIR}/../../data/config_nz.json COPYONLY)
configure_file(out.json ${CMAKE_CURRENT_BINARY_DIR}/../../data/out_nz.json COPYONLY)

# read the NewZealand grid and recharge
configure_file(data_5min/grid_nz.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/grid_nz.csv COPYONLY)
#configure_file(data_5min/recharge_nz.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/recharge_nz.csv COPYONLY)

# read the data for the model (nodes that are not in the grid read above are ignored)
configure_file(../../data_5min/lithology.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/lithology.csv COPYONLY)
configure_file(../../data_5min/elevation.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/elevation.csv COPYONLY)
configure_file(../../data_5min/river.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/river.csv COPYONLY)
configure_file(../../data_5min/initial_head.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/initial_head.csv COPYONLY)
configure_file(../../data_5min/surface_water_elevation.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/surface_water_elevation.csv COPYONLY)
configure_file(../../data_5min/recharge_withArcID.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/recharge_withArcID.csv COPYONLY)


## global and local wetlands
configure_file(../../data_5min/G_LOCWET.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/G_LOCWET.csv COPYONLY)
configure_file(../../data_5min/G_GLOWET.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/G_GLOWET.csv COPYONLY)
## global and local lakes
configure_file(../../data_5min/G_LOCLAK.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/G_LOCLAK.csv COPYONLY)
configure_file(../../data_5min/G_GLOLAK.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/G_GLOLAK.csv COPYONLY)
## slope
configure_file(../../data_5min/Slope/Africa_slope.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/Slope/Africa_slope.csv COPYONLY)
configure_file(../../data_5min/Slope/Australia_slope.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/Slope/Australia_slope.csv COPYONLY)
configure_file(../../data_5min/Slope/Eurasia_slope.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/Slope/Eurasia_slope.csv COPYONLY)
configure_file(../../data_5min/Slope/N_America_slope.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/Slope/N_America_slope.csv COPYONLY)
configure_file(../../data_5min/Slope/S_America_slope.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/Slope/S_America_slope.csv COPYONLY)
## e-folding
configure_file(../../data_5min/E-Folding/Africa_F.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/E-Folding/Africa_F.csv COPYONLY)
configure_file(../../data_5min/E-Folding/Australia_F.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/E-Folding/Australia_F.csv COPYONLY)
configure_file(../../data_5min/E-Folding/Eurasia_F.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/E-Folding/Eurasia_F.csv COPYONLY)
configure_file(../../data_5min/E-Folding/N_America_F.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/E-Folding/N_America_F.csv COPYONLY)
configure_file(../../data_5min/E-Folding/S_America_F.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/E-Folding/S_America_F.csv COPYONLY)
## mapping from spatID (used in most of the data) to arcID (used in recharge data)
configure_file(../../data_5min/spatID-arcID.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/spatID-arcID.csv COPYONLY)

