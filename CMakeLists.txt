## =================================================================== ##
## Copyright (C) 2014-2020 Kimura Ryo                                  ##
##                                                                     ##
## This Source Code Form is subject to the terms of the Mozilla Public ##
## License, v. 2.0. If a copy of the MPL was not distributed with this ##
## file, You can obtain one at http://mozilla.org/MPL/2.0/.            ##
## =================================================================== ##

cmake_minimum_required(VERSION 3.1.0)

project(libbsdf)

set_property(GLOBAL PROPERTY USE_FOLDERS ON)

set(HEADER_PATH include/libbsdf)
set(SOURCE_PATH src)

file(GLOB SOURCES
     "${HEADER_PATH}/*/*.h"
     "${SOURCE_PATH}/*/*.cpp"
     "${SOURCE_PATH}/*/*.c")

set(DIR_NAMES Brdf Common Reader Writer ReflectanceModel)
foreach(DIR_NAME ${DIR_NAMES})
    file(GLOB GROUP_HEADERS "${HEADER_PATH}/${DIR_NAME}/*.h")
    file(GLOB GROUP_SOURCES "${SOURCE_PATH}/${DIR_NAME}/*.cpp")
    source_group("${DIR_NAME}\\Header Files" FILES ${GROUP_HEADERS})
    source_group("${DIR_NAME}\\Source Files" FILES ${GROUP_SOURCES})
endforeach()

include_directories(include)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

# Enable all warnings
#add_definitions(-Wall)

if(EIGEN3_INCLUDE_DIR)
    include_directories("${EIGEN3_INCLUDE_DIR}")
endif()

if(MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP")
else()
    set (CMAKE_CXX_STANDARD 17)

    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -msse3")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse3")
endif()

option(USE_OpenMP "Use OpenMP" ON)
if(USE_OpenMP)
    find_package(OpenMP QUIET)
    if(OPENMP_FOUND)
        message("OpenMP found")
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    endif()
endif()

add_library(${PROJECT_NAME} ${SOURCES})

option(BUILD_LIBBSDF_APPLICATIONS "Enable to build libbsdf Applications" ON)
if(BUILD_LIBBSDF_APPLICATIONS)
    add_subdirectory(applications)

    if(NOT DEFINED LIBBSDF_CORE_FOLDER_NAME)
        set(LIBBSDF_CORE_FOLDER_NAME Core)
    endif()

    set_target_properties(${PROJECT_NAME} PROPERTIES FOLDER ${LIBBSDF_CORE_FOLDER_NAME})
endif()

target_link_libraries(${PROJECT_NAME})
