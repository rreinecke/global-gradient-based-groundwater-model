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
configure_file(config.json ${CMAKE_CURRENT_BINARY_DIR}/../../data/config.json COPYONLY)
configure_file(config_nz.json ${CMAKE_CURRENT_BINARY_DIR}/../../data/config_nz.json COPYONLY)
configure_file(config_gnc.json ${CMAKE_CURRENT_BINARY_DIR}/../../data/config_gnc.json COPYONLY)
configure_file(config_na.json ${CMAKE_CURRENT_BINARY_DIR}/../../data/config_na.json COPYONLY)
configure_file(out.json ${CMAKE_CURRENT_BINARY_DIR}/../../data/out.json COPYONLY)

# read the NewZealand grid and initial head
configure_file(data_5min/grid_nz.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/grid_nz.csv COPYONLY)
configure_file(data_5min/grid_nz_refined.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/grid_nz_refined.csv COPYONLY)
configure_file(data_5min/grid_nz_refined_16.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/grid_nz_refined_16.csv COPYONLY)
configure_file(data_5min/grid_na.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/grid_na.csv COPYONLY)
configure_file(../../data_5min/grid.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/grid.csv COPYONLY)
configure_file(../../data_5min/grid_refined_16.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/grid_refined_16.csv COPYONLY)
configure_file(../../data_5min/mean_wtd.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/mean_wtd.csv COPYONLY)

# read the data for the model (nodes that are not in the grid read above are ignored)
configure_file(../../data_5min/glhymps2.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/glhymps2.csv COPYONLY)
configure_file(../../data_5min/effective_porosity.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/effective_porosity.csv COPYONLY)
configure_file(../../data_5min/elevation.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/elevation.csv COPYONLY)
configure_file(../../data_5min/caspian_data.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/caspian_data.csv COPYONLY)

# rivers
configure_file(../../data_5min/river_extent.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/river_extent.csv COPYONLY)
configure_file(../../data_5min/surface_water_elevation.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/surface_water_elevation.csv COPYONLY)

# recharge
configure_file(../../data_5min/recharge_ISIMIP_2b.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/recharge_ISIMIP_2b.csv COPYONLY)

## global and local wetlands
configure_file(../../data_5min/G_LOCWET.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/G_LOCWET.csv COPYONLY)
configure_file(../../data_5min/G_GLOWET.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/G_GLOWET.csv COPYONLY)

## global and local lakes
configure_file(../../data_5min/G_LOCLAK.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/G_LOCLAK.csv COPYONLY)
configure_file(../../data_5min/G_GLOLAK.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/G_GLOLAK.csv COPYONLY)

## e-folding
configure_file(../../data_5min/E-Folding/Africa_F.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/E-Folding/Africa_F.csv COPYONLY)
configure_file(../../data_5min/E-Folding/Australia_F.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/E-Folding/Australia_F.csv COPYONLY)
configure_file(../../data_5min/E-Folding/Eurasia_F.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/E-Folding/Eurasia_F.csv COPYONLY)
configure_file(../../data_5min/E-Folding/N_America_F.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/E-Folding/N_America_F.csv COPYONLY)
configure_file(../../data_5min/E-Folding/S_America_F.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/E-Folding/S_America_F.csv COPYONLY)

## mapping from spatID (used in most of the data) to arcID (used in recharge data)
configure_file(../../data_5min/spatID-arcID.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/spatID-arcID.csv COPYONLY)

## zetas
configure_file(../../data_5min/Zetas/initial_zeta1_global.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/Zetas/initial_zeta1_global.csv COPYONLY)
configure_file(../../data_5min/Zetas/initial_zeta2_global.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/Zetas/initial_zeta2_global.csv COPYONLY)
configure_file(../../data_5min/Zetas/initial_zeta3_global.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/Zetas/initial_zeta3_global.csv COPYONLY)
configure_file(../../data_5min/Zetas/initial_zeta4_global.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/Zetas/initial_zeta4_global.csv COPYONLY)
