cmake_minimum_required(VERSION 3.15.3)
if(${CMAKE_VERSION} VERSION_LESS "3.15.3")
    message(STATUS "You need to upgrade to CMake 3.15.3")
endif()

project(GW)

# Build folder structure
add_subdirectory(src)
add_subdirectory(lib)

# CMAKE flags
set(CMAKE_MODULE_PATH "${GW_SOURCE_DIR}/cmake")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lpthread")

################# Libs ################
# Boost
set(BOOST_ROOT=/usr/local/lib)
set(Boost_NO_SYSTEM_PATHS=TRUE)
find_package(Boost 1.71.0 COMPONENTS system filesystem serialization log REQUIRED)
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
project(simple)
add_subdirectory(tests/SimpleModel)

set(simple_src ${SOURCE})
list(REMOVE_ITEM simple_src ${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp)

add_executable(simple_model tests/SimpleModel/simple.cpp ${simple_src} ${HEADERS})
target_link_libraries(simple_model pthread ${Boost_LIBRARIES} ${GMP_LIBRARIES})

################################
# Simple Demo Model for Variable Density Flow (VDF)
################################
project(simpleVDF)
add_subdirectory(tests/SimpleVDFModel)

set(simple_vdf_src ${SOURCE})
list(REMOVE_ITEM simple_vdf_src ${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp)
list(REMOVE_ITEM simple_vdf_src ${CMAKE_CURRENT_SOURCE_DIR}/tests/SimpleModel/simple.cpp)

add_executable(simple_vdf_model tests/SimpleVDFModel/simpleVDF.cpp ${simple_vdf_src} ${HEADERS})
target_link_libraries(simple_vdf_model pthread ${Boost_LIBRARIES} ${GMP_LIBRARIES})

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

# Prevent GoogleTest from overriding our compiler/linker options
# when building with Visual Studio
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

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

# Now simply link your own targets against gtest, gmock,
# etc. as appropriate

add_subdirectory(tests)
enable_testing()
#find_package(GTest REQUIRED)
#include_directories(${GTEST_INCLUDE_DIRS})
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DBOOST_LOG_DYN_LINK")

set(testing_hrd ${HEADERS})
set(testing_src ${SOURCE})
list(REMOVE_ITEM testing_src ${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp)
list(REMOVE_ITEM testing_src ${CMAKE_CURRENT_SOURCE_DIR}/tests/SimpleVDFModel/simpleVDF.cpp)


add_executable(runUnitTests ${testing_src} ${testing_hrd})
target_link_libraries(runUnitTests pthread ${Boost_LIBRARIES} ${GMP_LIBRARIES} gtest)
add_test(AllTests runUnitTests)



#project(Coverage)
#INCLUDE(cmake/CodeCoverage.cmake)
#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0 -fprofile-arcs -ftest-coverage")
#SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -O0 -fprofile-arcs -ftest-coverage")
#
#SETUP_TARGET_FOR_COVERAGE(
#		coverage_target     # Name for custom target.
#		runUnitTests        # Name of the test driver executable that runs the tests.
#		# NOTE! This should always have a ZERO as exit code
#		# otherwise the coverage generation will not complete.
#		coverage            # Name of output directory.
#)
