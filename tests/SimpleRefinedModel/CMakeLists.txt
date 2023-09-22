set(SOURCE
        ${SOURCE}
        ${CMAKE_CURRENT_SOURCE_DIR}/simpleRefined.cpp
        PARENT_SCOPE
        )
set(HEADERS
        ${HEADERS}
        ${CMAKE_CURRENT_SOURCE_DIR}/simpleRefined.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/SimpleRefinedDataReader.hpp
        PARENT_SCOPE
        )
# for simple model with refined nodes in grid center
configure_file(data/grid.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/grid_simpleRefined.csv COPYONLY)
configure_file(data/elevation.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/elevation_simpleRefined.csv COPYONLY)
configure_file(data/heads.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/heads_simpleRefined.csv COPYONLY)
configure_file(data/lithology.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/lithology_simpleRefined.csv COPYONLY)
configure_file(data/recharge.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/recharge_simpleRefined.csv COPYONLY)
configure_file(data/rivers.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/rivers_simpleRefined.csv COPYONLY)
configure_file(config.json ${CMAKE_CURRENT_BINARY_DIR}/../../data/config_simpleRefined.json COPYONLY)

# for simple model with refined nodes at coast
configure_file(data/grid_coast.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/grid_coast_simpleRefined.csv COPYONLY)
configure_file(data/elevation_coast.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/elevation_coast_simpleRefined.csv COPYONLY)
configure_file(data/heads_coast.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/heads_coast_simpleRefined.csv COPYONLY)
configure_file(data/lithology_coast.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/lithology_coast_simpleRefined.csv COPYONLY)
configure_file(data/recharge_coast.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/recharge_coast_simpleRefined.csv COPYONLY)
configure_file(data/rivers_coast.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/rivers_coast_simpleRefined.csv COPYONLY)
configure_file(data/ghb_coast.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/ghb_coast_simpleRefined.csv COPYONLY)
configure_file(config_coast.json ${CMAKE_CURRENT_BINARY_DIR}/../../data/config_coast_simpleRefined.json COPYONLY)

# output file
configure_file(out.json ${CMAKE_CURRENT_BINARY_DIR}/../../data/out_simpleRefined.json COPYONLY)