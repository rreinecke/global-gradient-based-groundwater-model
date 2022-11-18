cmake_minimum_required(VERSION 3.25.0)
if(${CMAKE_VERSION} VERSION_LESS "3.25.0")
    message(STATUS "You need to upgrade CMake! Minimum required version is 3.25.0.")
endif()

set(CMAKE_C_COMPILER "clang-13")
set(CMAKE_CXX_COMPILER "clang++-13")

# CMAKE flags
set(CMAKE_MODULE_PATH "${GW_SOURCE_DIR}/cmake")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++20")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99")

project(GW)

# Build folder structure
add_subdirectory(src)
add_subdirectory(lib)

set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lpthread")

################# Libs ################

# Boost
find_package(Boost 1.80.0 COMPONENTS system filesystem serialization log REQUIRED)
include_directories(${Boost_INCLUDE_DIR})
link_directories(${Boost_LIBRARY_DIR})
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DBOOST_LOG_DYN_LINK")

# OpenMP
FIND_PACKAGE(OpenMP)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -fopenmp")

# GMP
if (GMP_INCLUDE_DIR AND GMP_LIBRARIES)
    # Already in cache, be silent
    set(GMP_FIND_QUIETLY TRUE)
endif (GMP_INCLUDE_DIR AND GMP_LIBRARIES)
find_path(GMP_INCLUDE_DIR NAMES gmp.h)
find_library(GMP_LIBRARIES NAMES gmp libgmp)
find_library(GMPXX_LIBRARIES NAMES gmpxx libgmpxx)
MESSAGE(STATUS "GMP libs: " ${GMP_LIBRARIES} " " ${GMPXX_LIBRARIES})
include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(GMP DEFAULT_MSG GMP_INCLUDE_DIR GMP_LIBRARIES)

project(tests)
################################
# Simple Demo Model
################################
# project(simple)

# add_subdirectory(tests/SimpleModel)
 #add_executable(simple_model tests/SimpleModel/simple.cpp ${SOURCE} ${HEADERS})
# target_link_libraries(simple_model pthread ${Boost_LIBRARIES} ${GMP_LIBRARIES})

################################
# Simple Demo Model for Variable Density Flow (VDF)
################################
# project(simpleVDF)

# list(REMOVE_ITEM simple_src ${CMAKE_CURRENT_SOURCE_DIR}/tests/SimpleModel/simple.cpp)

#set(simple_vdf_hrd ${HEADERS})
#set(simple_vdf_src ${SOURCE})

#add_subdirectory(tests/SimpleVDFModel)
#add_executable(simple_vdf_model tests/SimpleVDFModel/simpleVDF.cpp ${SOURCE} ${HEADERS})
#target_link_libraries(simple_vdf_model pthread ${Boost_LIBRARIES} ${GMP_LIBRARIES})

################################
# New Zealand model
################################
# project(NewZealand)

# list(REMOVE_ITEM simple_vdf_model ${CMAKE_CURRENT_SOURCE_DIR}/tests/SimpleVDFModel/simpleVDF.cpp)

set(new_zealand_hrd ${HEADERS})
set(new_zealand_src ${SOURCE})

add_subdirectory(tests/NewZealand)
add_executable(new_zealand_model tests/NewZealand/runner.cpp ${SOURCE} ${HEADERS})
target_link_libraries(new_zealand_model pthread ${Boost_LIBRARIES} ${GMP_LIBRARIES})

################################
# Unit Tests
################################
project(Testing)

# Thanks to https://crascit.com/2015/07/25/cmake-gtest
# Download and unpack googletest at configure time
configure_file(CMakeLists_GTest.txt googletest-download/CMakeLists.txt)
execute_process(COMMAND "${CMAKE_COMMAND}" -G "${CMAKE_GENERATOR}" .
        WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/googletest-download" )
execute_process(COMMAND "${CMAKE_COMMAND}" --build .
        WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/googletest-download" )

# Add googletest directly to our build. This adds
# the following targets: gtest, gtest_main, gmock
# and gmock_main
add_subdirectory("${CMAKE_BINARY_DIR}/googletest-src"
        "${CMAKE_BINARY_DIR}/googletest-build")

# The gtest/gmock targets carry header search path
# dependencies automatically when using CMake 2.8.11 or
# later. Otherwise we have to add them here ourselves.
if(CMAKE_VERSION VERSION_LESS 2.8.11)
    include_directories("${gtest_SOURCE_DIR}/include"
            "${gmock_SOURCE_DIR}/include")
endif()

add_subdirectory(tests)
enable_testing()

set(testing_hrd ${HEADERS})
set(testing_src ${SOURCE})
list(REMOVE_ITEM testing_src ${CMAKE_CURRENT_SOURCE_DIR}/tests/SimpleModel/simple.cpp)
list(REMOVE_ITEM testing_src ${CMAKE_CURRENT_SOURCE_DIR}/tests/SimpleVDFModel/simpleVDF.cpp)
list(REMOVE_ITEM testing_src ${CMAKE_CURRENT_SOURCE_DIR}/tests/NewZealand/runner.cpp)


add_executable(runUnitTests ${testing_src} ${testing_hrd})
target_link_libraries(runUnitTests pthread ${Boost_LIBRARIES} ${GMP_LIBRARIES} gtest gtest_main gmock gmock_main)
add_test(AllTests runUnitTests)

################# Coverage ################
#project(Coverage)
#include(cmake/CodeCoverage.cmake)
#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0 -fprofile-arcs -ftest-coverage")
#SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -O0 -fprofile-arcs -ftest-coverage")
#
#set(COVERAGE_LCOV_EXCLUDES '${PROJECT_SOURCE_DIR}/lib/*' '/usr/include/*' '/usr/local/*' '*asio/*' '*googletest/*')

#SETUP_TARGET_FOR_COVERAGE_LCOV(
#                DEPENDENCIES runUnitTests
#                EXECUTABLE runUnitTests -j ${PROCESSOR_COUNT}
#                NAME coverage
#)