// =================================================================== //
// Copyright (C) 2015-2019 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#ifndef LIBBSDF_SDR_WRITER_H
#define LIBBSDF_SDR_WRITER_H

#include <string>

namespace lb {

class SampleSet2D;

/*!
 * \class SdrWriter
 * \brief The SdrWriter class provides the writer of a SDR and SDT file.
 */
class SdrWriter
{
public:
    /*! Writes specular reflectances or transmittances in a SDR or SDT file. */
    static bool write(const std::string& fileName,
                      const SampleSet2D& samples2D,
                      const std::string& comments = "");

    /*! Outputs character data of a SDR or SDT file to the stream. */
    static bool output(const SampleSet2D&   samples2D,
                       std::ostream&        stream,
                       const std::string&   comments = "");
};

} // namespace lb

#endif // LIBBSDF_SDR_WRITER_H
