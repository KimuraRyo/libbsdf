// =================================================================== //
// Copyright (C) 2015 Kimura Ryo                                       //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#ifndef LIBBSDF_FRESNEL_H
#define LIBBSDF_FRESNEL_H

#include <cmath>

namespace lb {

/*! Fresnel reflection. */
float fresnelReflection(float inTheta, float n1, float n2);

/*! Fresnel reflection. */
float fresnelReflection(float inTheta, float n);

/*
 * Implementation
 */

inline float fresnelReflection(float inTheta, float n1, float n2)
{
    float cosi = std::cos(inTheta);
    float sint = n1 / n2 * std::sin(inTheta);
    float cost = std::sqrt(1.0f - sint * sint);

    float rs = (n1 * cosi - n2 * cost) / (n1 * cosi + n2 * cost);
    float Rs = rs * rs;

    float rp = (n1 * cost - n2 * cosi) / (n1 * cost + n2 * cosi);
    float Rp = rp * rp;

    return (Rs + Rp) / 2.0f;
}

inline float fresnelReflection(float inTheta, float n)
{
    return fresnelReflection(inTheta, 1.0f, n);
}

} // namespace lb

#endif // LIBBSDF_FRESNEL_H
