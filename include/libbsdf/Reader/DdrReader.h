// =================================================================== //
// Copyright (C) 2014-2015 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#ifndef LIBBSDF_DDR_READER_H
#define LIBBSDF_DDR_READER_H

#include <string>

#include <libbsdf/Brdf/SpecularCoordinatesBrdf.h>

namespace lb {

/*!
 * \class DdrReader
 * \brief The DdrReader class provides the reader of a DDR and DDT file.
 *
 * A DDR (Diffuse Distribution Reflection) file includes BRDF, 
 * and a DDT (Diffuse Distribution Transparent) file includes BTDF.
 * It loads sample array and creates a lb::SpecularCoordinatesBrdf.
 */
class DdrReader
{
public:
    /*! Reads a DDR or DDT file and creates the BRDF of a specular coordinate system. */
    static SpecularCoordinatesBrdf* read(const std::string& fileName);
};

} // namespace lb

#endif // LIBBSDF_DDR_READER_H
