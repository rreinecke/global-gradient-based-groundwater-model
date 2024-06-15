set(SOURCE
        ${SOURCE}
        ${CMAKE_CURRENT_SOURCE_DIR}/simpleVDF.cpp
        PARENT_SCOPE
        )
set(HEADERS
        ${HEADERS}
        ${CMAKE_CURRENT_SOURCE_DIR}/simpleVDF.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/SimpleVDFDataReader.hpp
        PARENT_SCOPE
        )
configure_file(config.json ${CMAKE_CURRENT_BINARY_DIR}/../../data/config_simpleVDF.json COPYONLY)
configure_file(out.json ${CMAKE_CURRENT_BINARY_DIR}/../../data/out_simpleVDF.json COPYONLY)

# grid, elevation, initial zetas and initial head/wtd
configure_file(data/grid.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/grid_simpleVDF.csv COPYONLY) # done
configure_file(data/elevation.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/elevation_simpleVDF.csv COPYONLY)
configure_file(data/initial_zetas.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/initial_zetas_simpleVDF.csv COPYONLY)
configure_file(data/heads.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/heads_simpleVDF.csv COPYONLY)

# permeability, effective porosity, ghb, zones of sources and sinks
configure_file(data/ghb_conductance.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/ghb_conductance_simpleVDF.csv COPYONLY)
configure_file(data/zones_sources_sinks.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/zones_sources_sinks_simpleVDF.csv COPYONLY)

# recharge
configure_file(data/recharge.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/recharge_simpleVDF.csv COPYONLY)

# rivers
configure_file(data/rivers.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/rivers_simpleVDF.csv COPYONLY)