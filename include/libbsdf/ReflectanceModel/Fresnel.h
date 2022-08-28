// =================================================================== //
// Copyright (C) 2015-2021 Kimura Ryo                                  //
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
template <typename T>
T computeFresnel(const T& inTheta, const T& n1, const T& n2);

/*! Fresnel reflection. */
template <typename T>
T computeFresnel(const T& inTheta, const T& n);

/*! Fresnel reflection with a complex refractive index. */
template <typename T>
T computeComplexFresnel(const T& inTheta, const T& n, const T& k);

/*! Schlick's approximation of Fresnel reflection. */
template <typename T>
T computeSchlickFresnel(const T& inTheta, const T& n1, const T& n2 = T(1));

/*! Schlick's approximation of Fresnel reflection. */
template <typename ScalarT, typename ColorT>
ColorT computeSchlickFresnel(const ScalarT& dotNV, const ColorT& R0);

/*!
 * Spherical Gaussian approximation of Fresnel reflection.
 * Sébastien Lagarde, "Spherical Gaussian approximation for Blinn-Phong, Phong and Fresnel," June 2012.
 */
template <typename ScalarT, typename ColorT>
ColorT computeSphericalGaussianFresnel(const ScalarT& cosLN, const ColorT& R0);

/*
 * Implementation
 */

template <typename T>
T computeFresnel(const T& inTheta, const T& n1, const T& n2)
{
    assert(inTheta >= T(0) && inTheta <= T(PI_2_D));

    using std::cos;
    using std::sin;
    using std::sqrt;

    T cosI = cos(inTheta);
    T sinT = n1 / n2 * sin(inTheta);

    if (sinT >= T(1)) {
        // total internal reflection
        return T(1);
    }

    T cosT = sqrt(T(1) - sinT * sinT);

    T rs = (n1 * cosI - n2 * cosT) / (n1 * cosI + n2 * cosT);
    T Rs = rs * rs;

    T rp = (n1 * cosT - n2 * cosI) / (n1 * cosT + n2 * cosI);
    T Rp = rp * rp;

    return (Rs + Rp) / T(2);
}

template <typename T>
T computeFresnel(const T& inTheta, const T& n)
{
    return computeFresnel(inTheta, T(1), n);
}

template <typename T>
T computeComplexFresnel(const T& inTheta, const T& n, const T& k)
{
    assert(inTheta >= T(0) && inTheta <= T(PI_2_D));

    using std::cos;
    using std::sin;
    using std::sqrt;
    using std::tan;

    if (k == T(0)) {
        return computeFresnel(inTheta, n);
    }

    T sinI = sin(inTheta);
    T cosI = cos(inTheta);
    T tanI = tan(inTheta);

    T sqSinI = sinI * sinI;
    T sqCosI = cosI * cosI;
    T sqTanI = tanI * tanI;
    T sqN = n * n;
    T sqK = k * k;

    T nks = sqN - sqK - sqSinI;
    T sqA = (sqrt(nks * nks + T(4) * sqN * sqK) + (sqN - sqK - sqSinI)) / T(2);
    T sqB = (sqrt(nks * nks + T(4) * sqN * sqK) - (sqN - sqK - sqSinI)) / T(2);

    T a = sqrt(sqA);

    T Rs = (sqA + sqB - T(2) * a * cosI + sqCosI)
         / (sqA + sqB + T(2) * a * cosI + sqCosI);

    T Rp = Rs
         * (sqA + sqB - T(2) * a * sinI * tanI + sqSinI * sqTanI)
         / (sqA + sqB + T(2) * a * sinI * tanI + sqSinI * sqTanI);

    return (Rs + Rp) / T(2);
}

template <typename T>
T computeSchlickFresnel(const T& inTheta, const T& n1, const T& n2)
{
    assert(inTheta >= T(0) && inTheta <= T(PI_2_D));

    using std::cos;
    using std::pow;

    T r0 = (n1 - n2) / (n1 + n2);
    T R0 = r0 * r0;

    return R0 + (T(1) - R0) * pow(T(1) - cos(inTheta), T(5));
}

template <typename ScalarT, typename ColorT>
ColorT computeSchlickFresnel(const ScalarT& dotNV, const ColorT& R0)
{
    assert(dotNV >= ScalarT(0) && dotNV <= ScalarT(1));

    using std::pow;

    return R0 + (ColorT::Ones() - R0) * pow(ScalarT(1) - dotNV, ScalarT(5));
}

template <typename ScalarT, typename ColorT>
ColorT computeSphericalGaussianFresnel(const ScalarT& dotNV, const ColorT& R0)
{
    using std::exp2;

    return R0 + (ColorT::Ones() - R0) * exp2((ScalarT(-5.55473) * dotNV - ScalarT(6.98316)) * dotNV);
}

} // namespace lb

#endif // LIBBSDF_FRESNEL_H
