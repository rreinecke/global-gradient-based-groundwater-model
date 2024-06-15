set(SOURCE
        ${SOURCE}
        ${CMAKE_CURRENT_SOURCE_DIR}/usg_ex1.cpp
        PARENT_SCOPE
        )
set(HEADERS
        ${HEADERS}
        ${CMAKE_CURRENT_SOURCE_DIR}/usg_ex1.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/USG_ex1_DataReader.hpp
        PARENT_SCOPE
        )
# for simple model with refined nodes in grid center
configure_file(data/grid.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/grid_usg_ex1.csv COPYONLY)
configure_file(data/elevation.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/elevation_usg_ex1.csv COPYONLY)
configure_file(data/ghb_conductance.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/ghb_conductance_usg_ex1.csv COPYONLY)
configure_file(config.json ${CMAKE_CURRENT_BINARY_DIR}/../../data/config_usg_ex1.json COPYONLY)

# output file
configure_file(out.json ${CMAKE_CURRENT_BINARY_DIR}/../../data/out_usg_ex1.json COPYONLY)