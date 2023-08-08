#----------------------------------------------------------------------------
# Setup the project
#
cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
project(RadCorSub)

# Setup GNUInstallDirs for installation:
#----------------------------------------------------------------------------
include(GNUInstallDirs)
#----------------------------------------------------------------------------

#----------------------------------------------------------------------------
# test to link fortran function (check)

enable_language(Fortran)
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(dialect "-ffree-form -std=f2008 -fimplicit-none")
    set(bounds "-fbounds-check")
endif()
    set(CMAKE_Fortran_FLAGS "-ffree-line-length-none  -lowercase -extend_source -g -m64 -Wno-tabs -fno-f2c -DLinux -fno-automatic -finit-local-zero -ffixed-line-length-none -fno-second-underscore -funroll-loops")

#----------------------------------------------------------------------------
# The following lines are required so that the executable can dynamically locate the library dependencies independent of 
# use, i.e. don't skip the full RPATH for the build tree
SET(CMAKE_SKIP_BUILD_RPATH  FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_FULL_LIBDIR}")

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)


# the RPATH to be used when installing, but only if it's not a system directory
LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_FULL_LIBDIR}" isSystemDir)
IF("${isSystemDir}" STREQUAL "-1")
   SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_FULL_LIBDIR}")
ENDIF("${isSystemDir}" STREQUAL "-1")

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

#----------------------------------------------------------------------------

#
# Force use of c++11 standard (not sure)
#
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


#ROOT path
# commented out the following line so that this procedure will evolve with ROOT:
# set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/modules")
# You need to tell CMake where to find the ROOT installation. This can be done in a number of ways:
#   - ROOT built with classic configure/make use the provided $ROOTSYS/etc/cmake/FindROOT.cmake
#   - ROOT built with CMake. Add in CMAKE_PREFIX_PATH the installation prefix for ROOT
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED)

#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
#include(${ROOT_USE_FILE})

include_directories(${ROOT_INCLUDE_DIRS}
  ${PROJECT_SOURCE_DIR})

#include_directories(${PROJECT_SOURCE_DIR}/include)

message(PROJECT_SOURCE_DIR="${PROJECT_SOURCE_DIR}")

file(GLOB sources ${PROJECT_SOURCE_DIR}/*.f)

#${PROJECT_SOURCE_DIR}/src/*.cc)
#file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

add_executable(radcor_sbr radcor_wrap.cpp ${sources})
#${headers})

target_link_libraries(radcor_sbr ${ROOT_LIBRARIES} )


#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build N03. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
