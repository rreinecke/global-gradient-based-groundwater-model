cmake_minimum_required(VERSION 3.15.3)
add_subdirectory(Unit)
#add_subdirectory(Integration)
set(SOURCE
        ${SOURCE}
        ${CMAKE_CURRENT_SOURCE_DIR}/testrunner.cpp
        PARENT_SCOPE
        )
