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
configure_file(data/elevation.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/elevation_simpleVDF.csv COPYONLY)
configure_file(data/grid.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/grid_simpleVDF.csv COPYONLY)
configure_file(data/heads.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/heads_simpleVDF.csv COPYONLY)
configure_file(data/lithology.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/lithology_simpleVDF.csv COPYONLY)
configure_file(data/recharge.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/recharge_simpleVDF.csv COPYONLY)
configure_file(data/rivers.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/rivers_simpleVDF.csv COPYONLY)
configure_file(config.json ${CMAKE_CURRENT_BINARY_DIR}/../../data/config_simpleVDF.json COPYONLY)
configure_file(out.json ${CMAKE_CURRENT_BINARY_DIR}/../../data/out_simpleVDF.json COPYONLY)
