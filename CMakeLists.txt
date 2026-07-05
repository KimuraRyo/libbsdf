## =================================================================== ##
## Copyright (C) 2014-2026 Kimura Ryo                                  ##
##                                                                     ##
## This Source Code Form is subject to the terms of the Mozilla Public ##
## License, v. 2.0. If a copy of the MPL was not distributed with this ##
## file, You can obtain one at http://mozilla.org/MPL/2.0/.            ##
## =================================================================== ##

cmake_minimum_required(VERSION 3.22)

project(libbsdf VERSION 0.11.1)

set(CMAKE_CXX_STANDARD 17)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

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

# Set rgl-epfl/brdf-loader.
set(RGL_EPFL_BRDF_LOADER_PATH src/ThirdParty/rgl-epfl_brdf-loader)
set(RGL_EPFL_BRDF_LOADER_SOURCES
    "${RGL_EPFL_BRDF_LOADER_PATH}/powitacq.h"
    "${RGL_EPFL_BRDF_LOADER_PATH}/powitacq.inl"
    "${RGL_EPFL_BRDF_LOADER_PATH}/powitacq_rgb.h"
    "${RGL_EPFL_BRDF_LOADER_PATH}/powitacq_rgb.inl")
source_group("ThirdParty\\rgl-epfl_brdf-loader" FILES ${RGL_EPFL_BRDF_LOADER_SOURCES})
list(APPEND SOURCES ${RGL_EPFL_BRDF_LOADER_SOURCES})

set(CMAKE_INCLUDE_CURRENT_DIR ON)

find_package(Eigen3 REQUIRED)
if(Eigen3_FOUND)
    message("Found Eigen-${Eigen3_VERSION}: ${Eigen3_DIR}")
    set(EIGEN3_LIB "Eigen3::Eigen")
endif()

if(LIBBSDF_ENABLE_FITTING)
    find_package(Ceres REQUIRED)
    set(CERES_LIB "Ceres::ceres")
endif()

if(MSVC)
    add_compile_options(/W4)

    # Suppress warnings for rgl-epfl/brdf-loader headers
    set_source_files_properties(src/Reader/RglEpflBsdfReader.cpp PROPERTIES COMPILE_OPTIONS "/W2")

    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP")
else()
    if(NOT CMAKE_SYSTEM_PROCESSOR MATCHES "arm|aarch64|ARM64")
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -msse3")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse3")
    endif()
endif()

option(LIBBSDF_USE_OPENMP "Use OpenMP" ON)
if(LIBBSDF_USE_OPENMP)
    find_package(OpenMP QUIET)
    if(OPENMP_FOUND)
        message("Found OpenMP")
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    endif()
endif()

add_library(${PROJECT_NAME} ${SOURCES})
add_library(libbsdf::libbsdf ALIAS ${PROJECT_NAME})

target_include_directories(${PROJECT_NAME} PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

if(LIBBSDF_ENABLE_FITTING)
    target_compile_definitions(${PROJECT_NAME} PUBLIC LIBBSDF_USE_CERES_SOLVER)
endif()

option(LIBBSDF_BUILD_APPLICATIONS "Enable to build libbsdf Applications" ON)
if(LIBBSDF_BUILD_APPLICATIONS)
    add_subdirectory(applications)

    if(NOT DEFINED LIBBSDF_CORE_FOLDER_NAME)
        set(LIBBSDF_CORE_FOLDER_NAME Core)
    endif()

    set_target_properties(${PROJECT_NAME} PROPERTIES FOLDER ${LIBBSDF_CORE_FOLDER_NAME})
endif()

set(CMAKE_OPTIMIZE_DEPENDENCIES ON)
target_link_libraries(${PROJECT_NAME}
    PUBLIC
        Eigen3::Eigen
    PRIVATE
        ${CERES_LIB}
        $<$<BOOL:${OpenMP_CXX_FOUND}>:OpenMP::OpenMP_CXX>)

# ----------------------------------------------------------------------------
# Installation and package config for find_package(libbsdf).
# ----------------------------------------------------------------------------

install(TARGETS ${PROJECT_NAME}
    EXPORT libbsdfTargets
    ARCHIVE  DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY  DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME  DESTINATION ${CMAKE_INSTALL_BINDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(DIRECTORY include/libbsdf DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(EXPORT libbsdfTargets
    FILE      libbsdfTargets.cmake
    NAMESPACE libbsdf::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/libbsdf)

configure_package_config_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/libbsdfConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/libbsdfConfig.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/libbsdf)

write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/libbsdfConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion)

install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/libbsdfConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/libbsdfConfigVersion.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/libbsdf)
