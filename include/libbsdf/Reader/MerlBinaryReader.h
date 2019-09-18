// =================================================================== //
// Copyright (C) 2014-2015 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#ifndef LIBBSDF_MERL_BINARY_READER_H
#define LIBBSDF_MERL_BINARY_READER_H

#include <string>

#include <libbsdf/Brdf/HalfDifferenceCoordinatesBrdf.h>

namespace lb {

/*!
 * \class MerlBinaryReader
 * \brief The MerlBinaryReader class provides the reader for an isotropic BRDF measured by Matusik et al.
 *
 * lb::HalfDifferenceCoordinatesBrdf is created from loaded data.
 *
 * File format:
 * http://people.csail.mit.edu/wojciech/BRDFDatabase/code/
 *
 * Data:
 * http://www.merl.com/brdf/
 */
class MerlBinaryReader
{
public:
    /*! Reads a MERL binary file and creates the BRDF of a half difference coordinate system. */
    static HalfDifferenceCoordinatesBrdf* read(const std::string& fileName);
};

} // namespace lb

#endif // LIBBSDF_MERL_BINARY_READER_H
