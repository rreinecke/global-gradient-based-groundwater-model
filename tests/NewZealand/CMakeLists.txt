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
configure_file(data_5min/mother_NZ.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/mother_NZ_nz.csv COPYONLY)
configure_file(data_5min/litho_NZ.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/litho_NZ_nz.csv COPYONLY)
configure_file(data_5min/elevation_NZ.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/elevation_NZ_nz.csv COPYONLY)
configure_file(data_5min/recharge_NZ.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/recharge_NZ_nz.csv COPYONLY)
configure_file(data_5min/wtd_NZ.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/wtd_NZ_nz.csv COPYONLY)
configure_file(data_5min/blue_NZ.csv ${CMAKE_CURRENT_BINARY_DIR}/../../data/blue_NZ_nz.csv COPYONLY)
configure_file(out.json ${CMAKE_CURRENT_BINARY_DIR}/../../data/out_nz.json COPYONLY)
