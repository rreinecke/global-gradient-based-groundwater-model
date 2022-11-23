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
configure_file(data_5min/grid.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/grid_nz.csv COPYONLY)
configure_file(data_5min/lithology.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/lithology_nz.csv COPYONLY)
configure_file(data_5min/elevation.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/elevation_nz.csv COPYONLY)
configure_file(data_5min/recharge.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/recharge_nz.csv COPYONLY)
configure_file(data_5min/initial_head.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/initial_head_nz.csv COPYONLY)
configure_file(../data/G_LOCWET.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/G_LOCWET.csv COPYONLY)
configure_file(../data/G_GLOWET.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/G_GLOWET.csv COPYONLY)
configure_file(../data/G_LOCLAK.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/G_LOCLAK.csv COPYONLY)
configure_file(../data/G_GLOLAK.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/G_GLOLAK.csv COPYONLY)
configure_file(../data/Slope/Africa_slope.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/Slope/Africa_slope.csv COPYONLY)
configure_file(../data/Slope/Australia_slope.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/Slope/Australia_slope.csv COPYONLY)
configure_file(../data/Slope/Eurasia_slope.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/Slope/Eurasia_slope.csv COPYONLY)
configure_file(../data/Slope/N_America_slope.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/Slope/N_America_slope.csv COPYONLY)
configure_file(../data/Slope/S_America_slope.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/Slope/S_America_slope.csv COPYONLY)
configure_file(../data/E-Folding/Africa_F.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/E-Folding/Africa_F.csv COPYONLY)
configure_file(../data/E-Folding/Australia_F.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/E-Folding/Australia_F.csv COPYONLY)
configure_file(../data/E-Folding/Eurasia_F.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/E-Folding/Eurasia_F.csv COPYONLY)
configure_file(../data/E-Folding/N_America_F.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/E-Folding/N_America_F.csv COPYONLY)
configure_file(../data/E-Folding/S_America_F.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/E-Folding/S_America_F.csv COPYONLY)

configure_file(out.json ${CMAKE_CURRENT_BINARY_DIR}/../../data/out_nz.json COPYONLY)
