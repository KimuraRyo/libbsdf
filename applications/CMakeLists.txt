## =================================================================== ##
## Copyright (C) 2018-2019 Kimura Ryo                                  ##
##                                                                     ##
## This Source Code Form is subject to the terms of the Mozilla Public ##
## License, v. 2.0. If a copy of the MPL was not distributed with this ##
## file, You can obtain one at http://mozilla.org/MPL/2.0/.            ##
## =================================================================== ##

cmake_minimum_required(VERSION 3.1.0)

project(applications)

set(UTILITY_SOURCE_PATH Utility)
file(GLOB UTILITY_SOURCES
     "${UTILITY_SOURCE_PATH}/*.h"
     "${UTILITY_SOURCE_PATH}/*.cpp"
     "${UTILITY_SOURCE_PATH}/*.c")

if(NOT DEFINED LIBBSDF_APPLICATION_FOLDER_NAME)
    set(LIBBSDF_APPLICATION_FOLDER_NAME Applications)
endif()

include_directories(${UTILITY_SOURCE_PATH})
set(UTILITY_LIB_NAME AppUtility)
add_library(${UTILITY_LIB_NAME} ${UTILITY_SOURCES})
set_target_properties(${UTILITY_LIB_NAME} PROPERTIES FOLDER ${LIBBSDF_APPLICATION_FOLDER_NAME})

add_subdirectory(lbanalyze)
add_subdirectory(lbaverage)
add_subdirectory(lbcompare)
add_subdirectory(lbconv)
add_subdirectory(lbgen)
add_subdirectory(lbidentify)
add_subdirectory(lbinsert)
