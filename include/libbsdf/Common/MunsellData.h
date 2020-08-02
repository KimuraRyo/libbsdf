// =================================================================== //
// Copyright (C) 2020 Kimura Ryo                                       //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#ifndef LIBBSDF_MUNSELL_DATA_H
#define LIBBSDF_MUNSELL_DATA_H

namespace lb {

/*!
 * \struct MunsellData
 * \brief  The MunsellData struct provides Munsell renotation data.
 */
struct MunsellData
{
    static const char*  H[];    /*!< Hue array. */
    static const float  V[];    /*!< Value (lightness) array. */
    static const int    C[];    /*!< Chroma array. */
    static const float  xyY[];  /*!< xyY array. */

    static const int numColors = 4995; /*!< The number of colors. */

    /*!
     * The constant to normalize Y.
     * Y of white is more than 100, because magnesium oxide is used as the reference white.
     */
    static const float NORMALIZING_CONSTANT_Y;
};

} // namespace lb

#endif // LIBBSDF_MUNSELL_DATA_H
