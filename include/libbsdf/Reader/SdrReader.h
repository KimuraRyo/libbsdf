// =================================================================== //
// Copyright (C) 2014-2015 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#ifndef LIBBSDF_SDR_READER_H
#define LIBBSDF_SDR_READER_H

#include <string>

#include <libbsdf/Brdf/SampleSet2D.h>

namespace lb {

/*!
 * \class SdrReader
 * \brief The SdrReader class provides the reader of a SDR and SDT file.
 *
 * A SDR (Specular Distribution Reflection) file includes specular reflectances,
 * and a SDT (Specular Distribution Transparent) file includes specular transmittances.
 * It loads the list of reflectance and creates a lb::SampleSet2D.
 */
class SdrReader
{
public:
    /*! Reads a SDR or SDT file and creates sample points. */
    static SampleSet2D* read(const std::string& fileName);
};

} // namespace lb

#endif // LIBBSDF_SDR_READER_H
