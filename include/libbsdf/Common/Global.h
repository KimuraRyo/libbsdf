// =================================================================== //
// Copyright (C) 2014-2019 Kimura Ryo                                  //
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

namespace lb {

constexpr double PI_D   = 3.14159265358979323846;
constexpr double PI_2_D = 1.57079632679489661923;
constexpr double TAU_D  = 6.28318530717958647692;

constexpr float PI_F    = static_cast<float>(PI_D);
constexpr float PI_2_F  = static_cast<float>(PI_2_D);
constexpr float TAU_F   = static_cast<float>(TAU_D);

constexpr float EPSILON_F = std::numeric_limits<float>::epsilon();

/*! \brief The data type of a spectrum. */
using Spectrum = Eigen::ArrayXf;

/*! \brief The data type of spectra. */
using SpectrumList = std::vector<Spectrum, Eigen::aligned_allocator<Spectrum>>;

/*! \brief The output format of arrays and vectors. */
const Eigen::IOFormat LB_EIGEN_IO_FMT(Eigen::StreamPrecision, Eigen::DontAlignCols, " ", " ");

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
    GCMS4_FILE,
    INTEGRA_DDR_FILE,
    INTEGRA_DDT_FILE,
    INTEGRA_SDR_FILE,
    INTEGRA_SDT_FILE,
    LIGHTTOOLS_FILE,
    MERL_BINARY_FILE,
    ZEMAX_FILE
};

/*! \brief The data type of source. */
enum SourceType {
    UNKNOWN_SOURCE = 0,
    MEASURED_SOURCE,
    EDITED_SOURCE,
    GENERATED_SOURCE
};

} // namespace lb

#endif // LIBBSDF_GLOBAL_H
