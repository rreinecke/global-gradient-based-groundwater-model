set(SOURCE
        ${SOURCE}
        ${CMAKE_CURRENT_SOURCE_DIR}/simple.cpp
        PARENT_SCOPE
        )
set(HEADERS
        ${HEADERS}
        ${CMAKE_CURRENT_SOURCE_DIR}/simple.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/SimpleDataReader.hpp
        PARENT_SCOPE
        )
configure_file(data/elevation.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/elevation_simple.csv COPYONLY)
configure_file(data/grid.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/grid_simple.csv COPYONLY)
configure_file(data/heads.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/heads_simple.csv COPYONLY)
configure_file(data/lithology.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/lithology_simple.csv COPYONLY)
configure_file(data/recharge.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/recharge_simple.csv COPYONLY)
configure_file(data/rivers.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/rivers_simple.csv COPYONLY)
configure_file(config.json ${CMAKE_CURRENT_BINARY_DIR}/../../data/config_simple.json COPYONLY)
configure_file(out.json ${CMAKE_CURRENT_BINARY_DIR}/../../data/out_simple.json COPYONLY)