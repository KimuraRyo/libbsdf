// =================================================================== //
// Copyright (C) 2015-2024 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

/*!
 * \file    Version.h
 * \brief   The Version.h header file includes the version number of libbsdf.
 */

#ifndef LIBBSDF_VERSION_H
#define LIBBSDF_VERSION_H

#define LIBBSDF_MAJOR_VERSION 0
#define LIBBSDF_MINOR_VERSION 11
#define LIBBSDF_PATCH_VERSION 1

namespace lb {

/*! Gets the library version number. */
const char* getVersion();

} // namespace lb

#endif // LIBBSDF_VERSION_H
