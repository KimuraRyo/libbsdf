## =================================================================== ##
## Copyright (C) 2020-2022 Kimura Ryo                                  ##
##                                                                     ##
## This Source Code Form is subject to the terms of the Mozilla Public ##
## License, v. 2.0. If a copy of the MPL was not distributed with this ##
## file, You can obtain one at http://mozilla.org/MPL/2.0/.            ##
## =================================================================== ##

project(lbcombine)

set(SOURCES lbcombine.cpp)

include_directories(../${UTILITY_SOURCE_PATH})

add_executable(${PROJECT_NAME} ${SOURCES})
set_target_properties(${PROJECT_NAME} PROPERTIES FOLDER ${LIBBSDF_APPLICATION_FOLDER_NAME})

if(MSVC)
    set_target_properties(${PROJECT_NAME} PROPERTIES LINK_FLAGS "/SUBSYSTEM:CONSOLE")
endif()

target_link_libraries(${PROJECT_NAME} ${UTILITY_LIB_NAME} libbsdf)
