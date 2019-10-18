// =================================================================== //
// Copyright (C) 2014-2019 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#ifndef LIBBSDF_CIE_DATA_H
#define LIBBSDF_CIE_DATA_H

namespace lb {

/*!
 * \struct CieData
 * \brief  The CieData struct provides the data of CIE-XYZ, sRGB, and color matching functions.
 *
 * The range of wavelength is between 360nm and 830nm. The interval is 5nm.
 */
struct CieData
{
    static constexpr int numWavelengths = 95; /*!< The number of wavelengths: 95 samples. */

    static constexpr float minWavelength = 360.0f; /*!< The minimum value of wavelength: 360nm. */
    static constexpr float maxWavelength = 830.0f; /*!< The maximum value of wavelength: 830nm. */

    static const float XYZ[]; /*!< CIE(1931) 2-deg color matching function. */
    static const float D65[]; /*!< CIE Standard Illuminant D65 relative spectral power distribution. */

    /*! Transformation matrix from CIE-XYZ to sRGB. */
    static constexpr float XYZ_sRGB[9] =
    {
        3.2404542f, -1.5371385f, -0.4985314f,
        -0.9692660f, 1.8760108f, 0.0415560f,
        0.0556434f, -0.2040259f, 1.0572252f
    };

    /*! Transformation matrix from sRGB to CIE-XYZ. */
    static constexpr float sRGB_XYZ[9] =
    {
        0.4124564f, 0.3575761f, 0.1804375f,
        0.2126729f, 0.7151522f, 0.0721750f,
        0.0193339f, 0.1191920f, 0.9503041f
    };
};

} // namespace lb

#endif // LIBBSDF_CIE_DATA_H
