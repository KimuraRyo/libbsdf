// =================================================================== //
// Copyright (C) 2015-2017 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#ifndef LIBBSDF_FRESNEL_H
#define LIBBSDF_FRESNEL_H

#include <cmath>

#include <libbsdf/Common/Vector.h>

namespace lb {

/*! Fresnel reflection. */
float fresnel(float inTheta, float n1, float n2);

/*! Fresnel reflection. */
float fresnel(float inTheta, float n);

/*! Schlick's approximation. */
float schlickFresnel(float inTheta, float n1, float n2 = 1.0f);

/*! Schlick's approximation. */
Vec3 schlickFresnel(float cosLN, const Vec3& R0);

/*
 * Implementation
 */

inline float fresnel(float inTheta, float n1, float n2)
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

inline float fresnel(float inTheta, float n)
{
    return fresnel(inTheta, 1.0f, n);
}

inline float schlickFresnel(float inTheta, float n1, float n2)
{
    float r0 = (n1 - n2) / (n1 + n2);
    float R0 = r0 * r0;

    return R0 + (1.0f - R0) * std::pow(1.0f - std::cos(inTheta), 5.0f);
}

inline Vec3 schlickFresnel(float cosLN, const Vec3& R0)
{
    return R0 + (Vec3(1.0, 1.0, 1.0) - R0) * std::pow(1.0f - cosLN, 5.0f);
}

} // namespace lb

#endif // LIBBSDF_FRESNEL_H
