set(SOURCE
        ${SOURCE}
        ${CMAKE_CURRENT_SOURCE_DIR}/simpleVDF4.cpp
        PARENT_SCOPE
        )
set(HEADERS
        ${HEADERS}
        ${CMAKE_CURRENT_SOURCE_DIR}/simpleVDF4.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/SimpleVDF3DataReader.hpp
        PARENT_SCOPE
        )
configure_file(config.json ${CMAKE_CURRENT_BINARY_DIR}/../../data/config_simpleVDF3.json COPYONLY)
configure_file(out.json ${CMAKE_CURRENT_BINARY_DIR}/../../data/out_simpleVDF3.json COPYONLY)

# grid, elevation, initial zetas and initial head/wtd
configure_file(data/grid.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/grid_simpleVDF4.csv COPYONLY) # done
configure_file(data/elevation.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/elevation_simpleVDF4.csv COPYONLY)
configure_file(data/initial_zetas.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/initial_zetas_simpleVDF4.csv COPYONLY)
configure_file(data/heads.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/heads_simpleVDF4.csv COPYONLY)

# permeability and effective porosity
configure_file(data/unused/effective_porosity.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/effective_porosity_simpleVDF4.csv COPYONLY)
configure_file(data/lithology.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/lithology_simpleVDF4.csv COPYONLY)

# recharge
configure_file(data/recharge.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/recharge_simpleVDF4.csv COPYONLY)

# rivers
configure_file(data/rivers.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/rivers_simpleVDF4.csv COPYONLY)
#configure_file(data/surface_water_elevation.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/surface_water_elevation_simpleVDF4.csv COPYONLY)

## wetlands
configure_file(data/wetlands.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/wetlands_simpleVDF4.csv COPYONLY)
## lakes
configure_file(data/lakes.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/lakes_simpleVDF4.csv COPYONLY)
## e-folding
configure_file(data/E-Folding/efolding.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/E-Folding/efolding_simpleVDF4.csv COPYONLY)

## mapping from spatID (used in most of the data) to arcID (used in recharge data)
configure_file(data/spatID-arcID.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/spatID-arcID_simpleVDF4.csv COPYONLY)