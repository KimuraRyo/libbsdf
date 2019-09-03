// =================================================================== //
// Copyright (C) 2019 Kimura Ryo                                       //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#ifndef LIBBSDF_APPLICATIONS_UTILITY_H
#define LIBBSDF_APPLICATIONS_UTILITY_H

#include <iostream>
#include <string>

#include <libbsdf/Common/Utility.h>

namespace utility {

/*! Clamps the value of a parameter. */
template <typename T>
T clampParameter(const std::string& name, T val, T minVal, T maxVal);

/*
 * Implementation
 */

template <typename T>
T clampParameter(const std::string& name, T val, T minVal, T maxVal)
{
    if (val < minVal || val > maxVal) {
        val = clamp(val, minVal, maxVal);
        std::cout << "Warning: \"" + name + "\" is clamped to " << val << "." << std::endl;
    }

    return val;
}

} // namespace utility

#endif // LIBBSDF_APPLICATIONS_UTILITY_H
