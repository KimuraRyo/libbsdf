// =================================================================== //
// Copyright (C) 2014-2015 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

/*!
 * \file    DdrSdrUtility.h
 * \brief   Utility functions for DDR and DDT files.
 */

#ifndef LIBBSDF_DDR_SDR_UTILITY_H
#define LIBBSDF_DDR_SDR_UTILITY_H

#include <libbsdf/Reader/ReaderUtility.h>

namespace lb {
namespace ddr_sdr_utility {

enum SymmetryType {
    AXI_SYMMETRICAL,
    DIR_SYMMETRICAL,
    PLANE_SYMMETRICAL,
    ASYMMETRICAL,
    ASYMMETRICAL_4D
};

enum DataType {
    LUMINANCE_ABSOLUTE,
    INTENSITY_ABSOLUTE,
    LUMINANCE_RELATIVE,
    INTENSITY_RELATIVE
};

/*! Skips comment lines. */
void ignoreCommentLines(std::ifstream& fin);

} // namespace ddr_sdr_utility

inline void ddr_sdr_utility::ignoreCommentLines(std::ifstream& fin)
{
    reader_utility::ignoreCommentLines(fin, ";;");
}

} // namespace lb

#endif // LIBBSDF_DDR_SDR_UTILITY_H
