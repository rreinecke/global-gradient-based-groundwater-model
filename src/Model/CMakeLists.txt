set(SOURCE
        ${SOURCE}
        ${CMAKE_CURRENT_SOURCE_DIR}/Node.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/ExternalFlows.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/FluidMechanics.cpp
        PARENT_SCOPE
        )
set(HEADERS
        ${HEADERS}
        ${CMAKE_CURRENT_SOURCE_DIR}/PhysicalProperties.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/Node.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/Units.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/ExternalFlows.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/FluidMechanics.hpp
        PARENT_SCOPE
        )