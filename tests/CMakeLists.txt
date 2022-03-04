cmake_minimum_required(VERSION 2.8)
add_subdirectory(Unit)
add_subdirectory(files)
#add_subdirectory(Integration)
set(SOURCE
        ${SOURCE}
        ${CMAKE_CURRENT_SOURCE_DIR}/testrunner.cpp
        PARENT_SCOPE
        )
