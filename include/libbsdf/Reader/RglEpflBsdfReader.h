// =================================================================== //
// Copyright (C) 2026 Kimura Ryo                                       //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#ifndef LIBBSDF_RGL_EPFL_BSDF_READER_H
#define LIBBSDF_RGL_EPFL_BSDF_READER_H

#include <string>

#include <libbsdf/Brdf/SpecularCoordinatesBrdf.h>

namespace lb {

/*!
 * \class RglEpflBsdfReader
 * \brief The RglEpflBsdfReader class provides the reader for RGL-EPFL Material Database.
 *
 * lb::SpecularCoordinatesBrdf is created from loaded data.
 *
 * Data:
 * https://rgl.epfl.ch/materials
 */
class RglEpflBsdfReader
{
public:
    /*! Reads a RGL-EPFL BRDF file and creates the BRDF of a half difference coordinate system. */
    static SpecularCoordinatesBrdf* read(const std::string& fileName);
};

} // namespace lb

#endif // LIBBSDF_RGL_EPFL_BSDF_READER_H