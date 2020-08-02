// =================================================================== //
// Copyright (C) 2020 Kimura Ryo                                       //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#ifndef LIBBSDF_SSDD_UTILITY_H
#define LIBBSDF_SSDD_UTILITY_H

#include <Eigen/Core>

namespace lb {

namespace ssdd {
    constexpr char FILE_VERSION[] = "0.1";

    constexpr char VERSION[]    = "VERSION";
    constexpr char SOFTWARE[]   = "SOFTWARE";
    constexpr char API[]        = "API";
    constexpr char DATE[]       = "DATE";

    constexpr char DATA_TYPE[]                          = "DATA_TYPE";
    constexpr char DATA_TYPE_BRDF[]                     = "brdf";
    constexpr char DATA_TYPE_BTDF[]                     = "btdf";
    constexpr char DATA_TYPE_SPECULAR_REFLECTANCE[]     = "specular_reflectance";
    constexpr char DATA_TYPE_SPECULAR_TRANSMITTANCE[]   = "specular_transmittance";

    constexpr char COLOR_MODEL[]            = "COLOR_MODEL";
    constexpr char COLOR_MODEL_MONOCHROME[] = "monochrome";
    constexpr char COLOR_MODEL_RGB[]        = "rgb";
    constexpr char COLOR_MODEL_XYZ[]        = "xyz";
    constexpr char COLOR_MODEL_SPECTRUM[]   = "spectrum";

    constexpr char WAVELENGTH_LIST[] = "WAVELENGTH_LIST";

    constexpr char PARAM_TYPE[]             = "PARAM_TYPE";
    constexpr char PARAM_TYPE_HALF_DIFF[]   = "half_difference_coordinate_system";
    constexpr char PARAM_TYPE_SPECULAR[]    = "specular_coordinate_system";
    constexpr char PARAM_TYPE_SPHERICAL[]   = "spherical_coordinate_system";

    constexpr char REDUCTION_TYPE[]                     = "REDUCTION_TYPE";
    constexpr char REDUCTION_TYPE_BILATERAL_SYMMETRY[]  = "bilateral_symmetry";
    constexpr char REDUCTION_TYPE_RECIPROCITY[]         = "reciprocity";

    constexpr char PARAM0_LIST[] = "PARAM0_LIST";
    constexpr char PARAM1_LIST[] = "PARAM1_LIST";
    constexpr char PARAM2_LIST[] = "PARAM2_LIST";
    constexpr char PARAM3_LIST[] = "PARAM3_LIST";
    constexpr char PARAM4_LIST[] = "PARAM4_LIST";

    constexpr char SOURCE_TYPE[]            = "SOURCE_TYPE";
    constexpr char SOURCE_TYPE_EDITED[]     = "edited";
    constexpr char SOURCE_TYPE_GENERATED[]  = "generated";
    constexpr char SOURCE_TYPE_MEASURED[]   = "measured";

    constexpr char DEVICE[]             = "DEVICE";
    constexpr char CREATION_DATE[]      = "CREATION_DATE";
    constexpr char MEASUREMENT_DATE[]   = "MEASUREMENT_DATE";

    constexpr char DATA[]           = "DATA";
    constexpr char DATA_ASCII[]     = "ascii";
    constexpr char DATA_BINARY[]    = "binary";

    const Eigen::IOFormat LIST_FORMAT(Eigen::StreamPrecision, Eigen::DontAlignCols, " ", " ");
}

} // namespace lb

#endif // LIBBSDF_SSDD_UTILITY_H
