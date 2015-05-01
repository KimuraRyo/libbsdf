// =================================================================== //
// Copyright (C) 2015 Kimura Ryo                                       //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#ifndef LIBBSDF_ZEMAX_BSDF_READER_H
#define LIBBSDF_ZEMAX_BSDF_READER_H

#include <string>

#include <libbsdf/Brdf/SpecularCoordinatesBrdf.h>
#include <libbsdf/Reader/ReaderUtility.h>

namespace lb {

/*!
 * \class ZemaxBsdfReader
 * \brief The ZemaxBsdfReader class provides the reader for a Zemax BSDF file.
 *
 * A Zemax BSDF file can contain a BRDF or BTDF.
 *
 * File format:
 * https://www.zemax.com/support/knowledgebase/bsdf-data-interchange-file-format-specification
 */
class ZemaxBsdfReader
{
public:
    /*!
     * Reads a Zemax BSDF file and creates the two-sided material of a spherical coordinate system.
     * \a dataType returns lb::BRDF_DATA or lb::BTDF_DATA.
     */
    static SpecularCoordinatesBrdf* read(const std::string& fileName, DataType* dataType);

private:
    enum SymmetryType {
        UNKNOWN_SYMMETRY = 0,
        PLANE_SYMMETRICAL,
        ASYMMETRICAL,
        ASYMMETRICAL_4D
    };

    /*! Skips comment lines. */
    static void ignoreCommentLines(std::ifstream& fin);
};

inline void ZemaxBsdfReader::ignoreCommentLines(std::ifstream& fin)
{
    reader_utility::ignoreCommentLines(fin, "#");
}

} // namespace lb

#endif // LIBBSDF_ZEMAX_BSDF_READER_H
