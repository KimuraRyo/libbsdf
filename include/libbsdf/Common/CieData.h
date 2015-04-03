// =================================================================== //
// Copyright (C) 2014-2015 Kimura Ryo                                  //
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
    static const int numWavelengths; /*!< The number of wavelengths: 95 samples. */
    
    static const float minWavelength; /*!< The minimum value of wavelength: 360nm. */
    static const float maxWavelength; /*!< The maximum value of wavelength: 830nm. */

    static const float XYZ[]; /*!< CIE(1931) 2-deg color matching function. */
    static const float D65[]; /*!< CIE Standard Illuminant D65 relative spectral power distribution. */

    static const float XYZ_sRGB[9]; /*!< Transformation matrix from CIE-XYZ to sRGB. */
};

} // namespace lb

#endif // LIBBSDF_CIE_DATA_H
