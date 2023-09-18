set(SOURCE
        ${SOURCE}
        ${CMAKE_CURRENT_SOURCE_DIR}/swi2_ex2.cpp
        PARENT_SCOPE
        )
set(HEADERS
        ${HEADERS}
        ${CMAKE_CURRENT_SOURCE_DIR}/swi2_ex2.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/SWI2_ex2_DataReader.hpp
        PARENT_SCOPE
        )
configure_file(data/ghb_conductance.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/ghb_conductance_swi2_ex2.csv COPYONLY)
configure_file(data/elevation.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/elevation_swi2_ex2.csv COPYONLY)
configure_file(data/grid.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/grid_swi2_ex2.csv COPYONLY)
configure_file(data/unused/initial_head.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/initial_head_swi2_ex2.csv COPYONLY)
configure_file(data/initial_zetas.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/Zetas/initial_zetas_swi2_ex2.csv COPYONLY)
configure_file(data/effective_porosity.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/effective_porosity_swi2_ex2.csv COPYONLY)
configure_file(data/zones_sources_sinks.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/zones_sources_sinks_swi2_ex2.csv COPYONLY)
configure_file(data/unused/lithology.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/lithology_swi2_ex2.csv COPYONLY)
configure_file(data/unused/recharge.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/recharge_swi2_ex2.csv COPYONLY)
configure_file(config.json ${CMAKE_CURRENT_BINARY_DIR}/../../data/config_swi2_ex2.json COPYONLY)
configure_file(out.json ${CMAKE_CURRENT_BINARY_DIR}/../../data/out_swi2_ex2.json COPYONLY)