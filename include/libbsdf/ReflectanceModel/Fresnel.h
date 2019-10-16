// =================================================================== //
// Copyright (C) 2015-2019 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#ifndef LIBBSDF_FRESNEL_H
#define LIBBSDF_FRESNEL_H

#define _USE_MATH_DEFINES
#include <cmath>

#include <libbsdf/Common/Utility.h>

namespace lb {

/*! Fresnel reflection. */
float fresnel(float inTheta, float n1, float n2);

/*! Fresnel reflection. */
float fresnel(float inTheta, float n);

/*! Fresnel reflection with a complex refractive index. */
float fresnelComplex(float inTheta, float n, float k);

/*! Schlick's approximation of Fresnel reflection. */
float fresnelSchlick(float inTheta, float n1, float n2 = 1.0f);

/*! Schlick's approximation of Fresnel reflection. */
Vec3 fresnelSchlick(float cosLN, const Vec3& R0);

/*
 * Implementation
 */

inline float fresnel(float inTheta, float n1, float n2)
{
    assert(inTheta >= 0.0f && inTheta <= PI_2_F);

    float cosI = std::cos(inTheta);
    float sinT = n1 / n2 * std::sin(inTheta);

    if (sinT >= 1.0f) {
        // total internal reflection
        return 1.0f;
    }

    float cosT = std::sqrt(1.0f - sinT * sinT);

    float rs = (n1 * cosI - n2 * cosT) / (n1 * cosI + n2 * cosT);
    float Rs = rs * rs;

    float rp = (n1 * cosT - n2 * cosI) / (n1 * cosT + n2 * cosI);
    float Rp = rp * rp;

    return (Rs + Rp) / 2.0f;
}

inline float fresnel(float inTheta, float n)
{
    return fresnel(inTheta, 1.0f, n);
}

inline float fresnelComplex(float inTheta, float n, float k)
{
    if (k == 0.0f) {
        return fresnel(inTheta, n);
    }

    assert(inTheta >= 0.0f && inTheta <= PI_2_F);

    using std::sqrt;

    float sinI = std::sin(inTheta);
    float cosI = std::cos(inTheta);
    float tanI = std::tan(inTheta);

    float sqSinI = sinI * sinI;
    float sqCosI = cosI * cosI;
    float sqTanI = tanI * tanI;
    float sqN = n * n;
    float sqK = k * k;

    float nks = sqN - sqK - sqSinI;
    float sqA = (sqrt(nks * nks + 4.0f * sqN * sqK) + (sqN - sqK - sqSinI)) / 2.0f;
    float sqB = (sqrt(nks * nks + 4.0f * sqN * sqK) - (sqN - sqK - sqSinI)) / 2.0f;

    float a = sqrt(sqA);

    float Rs = (sqA + sqB - 2.0f * a * cosI + sqCosI)
             / (sqA + sqB + 2.0f * a * cosI + sqCosI);

    float Rp = Rs
             * (sqA + sqB - 2.0f * a * sinI * tanI + sqSinI * sqTanI)
             / (sqA + sqB + 2.0f * a * sinI * tanI + sqSinI * sqTanI);

    return (Rs + Rp) / 2.0f;
}

inline float fresnelSchlick(float inTheta, float n1, float n2)
{
    assert(inTheta >= 0.0f && inTheta <= PI_2_F);

    float r0 = (n1 - n2) / (n1 + n2);
    float R0 = r0 * r0;

    return R0 + (1.0f - R0) * std::pow(1.0f - std::cos(inTheta), 5.0f);
}

inline Vec3 fresnelSchlick(float cosLN, const Vec3& R0)
{
    return R0 + (Vec3(1.0, 1.0, 1.0) - R0) * std::pow(1.0f - cosLN, 5.0f);
}

} // namespace lb

#endif // LIBBSDF_FRESNEL_H
