## =================================================================== ##
## Copyright (C) 2014-2022 Kimura Ryo                                  ##
##                                                                     ##
## This Source Code Form is subject to the terms of the Mozilla Public ##
## License, v. 2.0. If a copy of the MPL was not distributed with this ##
## file, You can obtain one at http://mozilla.org/MPL/2.0/.            ##
## =================================================================== ##

cmake_minimum_required(VERSION 3.5)

project(libbsdf)

set_property(GLOBAL PROPERTY USE_FOLDERS ON)

option(LIBBSDF_ENABLE_FITTING "Enable fitting functions with Ceres Solver" OFF)

set(HEADER_PATH include/libbsdf)
set(SOURCE_PATH src)

# Create groups for each directory and add file paths to the list of source code.
macro(set_sources DIR_NAME)
    file(GLOB GROUP_HEADERS "${HEADER_PATH}/${DIR_NAME}/*.h")
    file(GLOB GROUP_SOURCES "${SOURCE_PATH}/${DIR_NAME}/*.cpp")
    source_group("${DIR_NAME}\\Header Files" FILES ${GROUP_HEADERS})
    source_group("${DIR_NAME}\\Source Files" FILES ${GROUP_SOURCES})

    list(APPEND SOURCES ${GROUP_HEADERS} ${GROUP_SOURCES})
endmacro()

set(SOURCES)
set_sources(Brdf)
set_sources(Common)
if(LIBBSDF_ENABLE_FITTING)
    set_sources(Fitter)
endif()
set_sources(Reader)
set_sources(ReflectanceModel)
set_sources(Writer)

# Set delaunator-cpp.
set(DELAUNATOR_PATH src/ThirdParty/delaunator-cpp)
set(DELAUNATOR_SOURCES "${DELAUNATOR_PATH}/delaunator.hpp" "${DELAUNATOR_PATH}/delaunator.cpp")
source_group("ThirdParty\\delaunator-cpp" FILES ${DELAUNATOR_SOURCES})
list(APPEND SOURCES ${DELAUNATOR_SOURCES})

include_directories(include)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

if(EIGEN3_INCLUDE_DIR)
    include_directories("${EIGEN3_INCLUDE_DIR}")
endif()

if(LIBBSDF_ENABLE_FITTING)
    find_package(Ceres REQUIRED)
    set(USED_LIBRARIES "Ceres::ceres")
    add_definitions(-DLIBBSDF_USE_CERES_SOLVER)
endif()

if(MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP")
else()
    set (CMAKE_CXX_STANDARD 17)

    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -msse3")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse3")
endif()

option(LIBBSDF_USE_OPENMP "Use OpenMP" ON)
if(LIBBSDF_USE_OPENMP)
    find_package(OpenMP QUIET)
    if(OPENMP_FOUND)
        message("OpenMP found")
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    endif()
endif()

add_library(${PROJECT_NAME} ${SOURCES})

option(LIBBSDF_BUILD_APPLICATIONS "Enable to build libbsdf Applications" ON)
if(LIBBSDF_BUILD_APPLICATIONS)
    add_subdirectory(applications)

    if(NOT DEFINED LIBBSDF_CORE_FOLDER_NAME)
        set(LIBBSDF_CORE_FOLDER_NAME Core)
    endif()

    set_target_properties(${PROJECT_NAME} PROPERTIES FOLDER ${LIBBSDF_CORE_FOLDER_NAME})
endif()

target_link_libraries(${PROJECT_NAME} ${USED_LIBRARIES})
