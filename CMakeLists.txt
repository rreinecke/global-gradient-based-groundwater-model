cmake_minimum_required(VERSION 3.25.0)
if(${CMAKE_VERSION} VERSION_LESS "3.25.0")
    message(STATUS "You need to upgrade CMake! Minimum required version is 3.25.0.")
endif()

set(CMAKE_C_COMPILER "clang-14")
set(CMAKE_CXX_COMPILER "clang++-14")

# CMAKE flags
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++20")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99")

project(3GM)

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

################################
# Simple Demo Model
################################
add_subdirectory(tests/SimpleModel)
add_executable(simple_model tests/SimpleModel/simple.cpp ${SOURCE} ${HEADERS})
target_link_libraries(simple_model pthread ${Boost_LIBRARIES} ${GMP_LIBRARIES})

################################
# Simple Demo Models for refined grids
################################
add_subdirectory(tests/USG_ex1)
set(usg_ex1_src ${SOURCE})
list(REMOVE_ITEM usg_ex1_src ${CMAKE_SOURCE_DIR}/tests/SimpleModel/simple.cpp)
add_executable(usg_ex1 tests/USG_ex1/usg_ex1.cpp ${usg_ex1_src} ${HEADERS})
target_link_libraries(usg_ex1 pthread ${Boost_LIBRARIES} ${GMP_LIBRARIES})

add_subdirectory(tests/SimpleRefinedModel)
set(simple_refined_src ${SOURCE})
list(REMOVE_ITEM simple_refined_src ${CMAKE_SOURCE_DIR}/tests/SimpleModel/simple.cpp)
list(REMOVE_ITEM simple_refined_src ${CMAKE_SOURCE_DIR}/tests/USG_ex1/usg_ex1.cpp)
add_executable(simple_refined_model tests/SimpleRefinedModel/simpleRefined.cpp ${simple_refined_src} ${HEADERS})
target_link_libraries(simple_refined_model pthread ${Boost_LIBRARIES} ${GMP_LIBRARIES})

################################
# Simple Demo Models for Variable Density Flow (VDF)
################################
add_subdirectory(tests/SWI2_ex1)
set(swi2_ex1_src ${SOURCE})
list(REMOVE_ITEM swi2_ex1_src ${CMAKE_CURRENT_SOURCE_DIR}/tests/SimpleModel/simple.cpp)
list(REMOVE_ITEM swi2_ex1_src ${CMAKE_SOURCE_DIR}/tests/USG_ex1/usg_ex1.cpp)
list(REMOVE_ITEM swi2_ex1_src ${CMAKE_CURRENT_SOURCE_DIR}/tests/SimpleRefinedModel/simpleRefined.cpp)
add_executable(swi2_ex1 tests/SWI2_ex1/swi2_ex1.cpp ${swi2_ex1_src} ${HEADERS})
target_link_libraries(swi2_ex1 pthread ${Boost_LIBRARIES} ${GMP_LIBRARIES})

add_subdirectory(tests/SWI2_ex2)
set(swi2_ex2_src ${SOURCE})
list(REMOVE_ITEM swi2_ex2_src ${CMAKE_SOURCE_DIR}/tests/SimpleModel/simple.cpp)
list(REMOVE_ITEM swi2_ex2_src ${CMAKE_SOURCE_DIR}/tests/USG_ex1/usg_ex1.cpp)
list(REMOVE_ITEM swi2_ex2_src ${CMAKE_SOURCE_DIR}/tests/SimpleRefinedModel/simpleRefined.cpp)
list(REMOVE_ITEM swi2_ex2_src ${CMAKE_SOURCE_DIR}/tests/SWI2_ex1/swi2_ex1.cpp)
add_executable(swi2_ex2 tests/SWI2_ex2/swi2_ex2.cpp ${swi2_ex2_src} ${HEADERS})
target_link_libraries(swi2_ex2 pthread ${Boost_LIBRARIES} ${GMP_LIBRARIES})

add_subdirectory(tests/SWI2_ex3)
set(swi2_ex3_src ${SOURCE})
list(REMOVE_ITEM swi2_ex3_src ${CMAKE_SOURCE_DIR}/tests/SimpleModel/simple.cpp)
list(REMOVE_ITEM swi2_ex3_src ${CMAKE_SOURCE_DIR}/tests/USG_ex1/usg_ex1.cpp)
list(REMOVE_ITEM swi2_ex3_src ${CMAKE_SOURCE_DIR}/tests/SimpleRefinedModel/simpleRefined.cpp)
list(REMOVE_ITEM swi2_ex3_src ${CMAKE_SOURCE_DIR}/tests/SWI2_ex1/swi2_ex1.cpp)
list(REMOVE_ITEM swi2_ex3_src ${CMAKE_SOURCE_DIR}/tests/SWI2_ex2/swi2_ex2.cpp)
add_executable(swi2_ex3 tests/SWI2_ex3/swi2_ex3.cpp ${swi2_ex3_src} ${HEADERS})
target_link_libraries(swi2_ex3 pthread ${Boost_LIBRARIES} ${GMP_LIBRARIES})

add_subdirectory(tests/SimpleVDFModel)
set(simple_vdf_src ${SOURCE})
list(REMOVE_ITEM simple_vdf_src ${CMAKE_SOURCE_DIR}/tests/SimpleModel/simple.cpp)
list(REMOVE_ITEM simple_vdf_src ${CMAKE_SOURCE_DIR}/tests/USG_ex1/usg_ex1.cpp)
list(REMOVE_ITEM simple_vdf_src ${CMAKE_SOURCE_DIR}/tests/SimpleRefinedModel/simpleRefined.cpp)
list(REMOVE_ITEM simple_vdf_src ${CMAKE_SOURCE_DIR}/tests/SWI2_ex1/swi2_ex1.cpp)
list(REMOVE_ITEM simple_vdf_src ${CMAKE_SOURCE_DIR}/tests/SWI2_ex2/swi2_ex2.cpp)
list(REMOVE_ITEM simple_vdf_src ${CMAKE_SOURCE_DIR}/tests/SWI2_ex3/swi2_ex3.cpp)
add_executable(simple_vdf_model tests/SimpleVDFModel/simpleVDF.cpp ${simple_vdf_src} ${HEADERS})
target_link_libraries(simple_vdf_model pthread ${Boost_LIBRARIES} ${GMP_LIBRARIES})

################################
# Large scale model (with Variable Density Flow (VDF))
################################
add_subdirectory(tests/LargeModel)
set(large_model_src ${SOURCE})
list(REMOVE_ITEM large_model_src ${CMAKE_SOURCE_DIR}/tests/SimpleModel/simple.cpp)
list(REMOVE_ITEM large_model_src ${CMAKE_SOURCE_DIR}/tests/USG_ex1/usg_ex1.cpp)
list(REMOVE_ITEM large_model_src ${CMAKE_SOURCE_DIR}/tests/SimpleRefinedModel/simpleRefined.cpp)
list(REMOVE_ITEM large_model_src ${CMAKE_SOURCE_DIR}/tests/SWI2_ex1/swi2_ex1.cpp)
list(REMOVE_ITEM large_model_src ${CMAKE_SOURCE_DIR}/tests/SWI2_ex2/swi2_ex2.cpp)
list(REMOVE_ITEM large_model_src ${CMAKE_SOURCE_DIR}/tests/SWI2_ex3/swi2_ex3.cpp)
list(REMOVE_ITEM large_model_src ${CMAKE_SOURCE_DIR}/tests/SimpleVDFModel/simpleVDF.cpp)
add_executable(large_model tests/LargeModel/runner.cpp ${large_model_src} ${HEADERS})
target_link_libraries(large_model pthread ${Boost_LIBRARIES} ${GMP_LIBRARIES})

################################
# Unit Tests
################################
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
list(REMOVE_ITEM testing_src ${CMAKE_SOURCE_DIR}/tests/SimpleModel/simple.cpp)
list(REMOVE_ITEM testing_src ${CMAKE_SOURCE_DIR}/tests/USG_ex1/usg_ex1.cpp)
list(REMOVE_ITEM testing_src ${CMAKE_SOURCE_DIR}/tests/SimpleRefinedModel/simpleRefined.cpp)
list(REMOVE_ITEM testing_src ${CMAKE_SOURCE_DIR}/tests/SWI2_ex1/swi2_ex1.cpp)
list(REMOVE_ITEM testing_src ${CMAKE_SOURCE_DIR}/tests/SWI2_ex2/swi2_ex2.cpp)
list(REMOVE_ITEM testing_src ${CMAKE_SOURCE_DIR}/tests/SWI2_ex3/swi2_ex3.cpp)
list(REMOVE_ITEM testing_src ${CMAKE_SOURCE_DIR}/tests/SimpleVDFModel/simpleVDF.cpp)
list(REMOVE_ITEM testing_src ${CMAKE_SOURCE_DIR}/tests/LargeModel/runner.cpp)

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