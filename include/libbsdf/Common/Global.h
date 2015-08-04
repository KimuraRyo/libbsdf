// =================================================================== //
// Copyright (C) 2014-2015 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

/*!
 * \file    Global.h
 * \brief   The Global.h header file includes the global declarations.
 */

#ifndef LIBBSDF_GLOBAL_H
#define LIBBSDF_GLOBAL_H

#include <vector>

#include <Eigen/Core>

#include <libbsdf/Common/Array.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace lb {

const float PI_F = 3.14159265358979323846f;
const float PI_2_F = 1.57079632679489661923f;
const float EPSILON_F = std::numeric_limits<float>::epsilon();

/*! \brief The data type of a spectrum. */
typedef Arrayf Spectrum;

/*! \brief The data type of spectra. */
typedef std::vector<Spectrum, Eigen::aligned_allocator<Spectrum> > SpectrumList;

/*! \brief The color models and spaces of the data in spectra. */
enum ColorModel {
    UNKNOWN_MODEL = 0,
    MONOCHROMATIC_MODEL,
    RGB_MODEL,
    XYZ_MODEL,
    SPECTRAL_MODEL
};

/*! \brief The data type of scatter. */
enum DataType {
    UNKNOWN_DATA = 0,
    BRDF_DATA,
    BTDF_DATA,
    SPECULAR_REFLECTANCE_DATA,
    SPECULAR_TRANSMITTANCE_DATA
};

/*! \brief The type of data files. */
enum FileType {
    UNKNOWN_FILE = 0,
    ASTM_FILE,
    INTEGRA_DDR_FILE,
    INTEGRA_DDT_FILE,
    INTEGRA_SDR_FILE,
    INTEGRA_SDT_FILE,
    LIGHTTOOLS_FILE,
    MERL_BINARY_FILE,
    ZEMAX_FILE
};

} // namespace lb

#endif // LIBBSDF_GLOBAL_H
