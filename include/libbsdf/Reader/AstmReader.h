// =================================================================== //
// Copyright (C) 2014-2015 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#ifndef LIBBSDF_ASTM_READER_H
#define LIBBSDF_ASTM_READER_H

#include <string>

#include <libbsdf/Brdf/SphericalCoordinatesBrdf.h>

namespace lb {

/*!
 * \class AstmReader
 * \brief The AstmReader class provides the reader for an ASTM E1392-96(2002) file.
 *
 * lb::SphericalCoordinatesBrdf is created from loaded data.
 * 
 * File format:
 * http://www.astm.org/Standards/E1392.htm
 */
class AstmReader
{
public:
    /*! Reads an ASTM file and creates the BRDF of a spherical coordinate system. */
    static SphericalCoordinatesBrdf* read(const std::string& fileName);
};

} // namespace lb

#endif // LIBBSDF_ASTM_READER_H
