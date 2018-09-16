// =================================================================== //
// Copyright (C) 2014-2018 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

/*!
 * \file    Utility.h
 * \brief   The Utility.h header file includes the utility functions for libbsdf.
 */

#ifndef LIBBSDF_UTILITY_H
#define LIBBSDF_UTILITY_H

#include <iostream>

#include <libbsdf/Common/CentripetalCatmullRomSpline.h>
#include <libbsdf/Common/CieData.h>
#include <libbsdf/Common/Global.h>
#include <libbsdf/Common/Vector.h>

namespace lb {

/*! \brief Clamps a value between a minimum and maximum value. */
template <typename T>
T clamp(T value, T minValue, T maxValue);

/*! \brief Returns true if two values are nearly equal. */
template <typename T>
bool isEqual(T lhs, T rhs);

/*! \brief Returns -1 for a negative number, +1 for a positive number, and 0 for 0. */
template <typename T> int sign(T val) {
    return (T(0) < val) - (val < T(0));
}

/*! \brief Computes linearly-interpolated values. */
template <typename T>
T lerp(const T& v0, const T& v1, float t);

/*! \brief Interpolates smoothly between two input values with cubic Hermite interpolation. */
template <typename T>
T smoothstep(const T& v0, const T& v1, const T& t);

/*! \brief Interpolates smoothly between two input values with 5th-order Hermite interpolation. */
template <typename T>
T smootherstep(const T& v0, const T& v1, const T& t);

/*! \brief Computes smoothly interpolated values with cubic Hermite interpolation. */
template <typename T>
T hermiteInterpolation3(const T& v0, const T& v1, float t);

/*! \brief Computes smoothly interpolated values with 5th-order Hermite interpolation. */
template <typename T>
T hermiteInterpolation5(const T& v0, const T& v1, float t);

/*! \brief Computes interpolated values using uniform Catmull-Rom spline. */
template <typename T>
T catmullRomSpline(const T& v0, const T& v1, const T& v2, const T& v3, float t);

/*!
 * \brief Computes a interpolated value using centripetal Catmull-Rom spline
 *        at \a pos in [\a pos1,\a pos2].
 */
template <typename T>
T catmullRomSpline(const T& pos0, const T& pos1, const T& pos2, const T& pos3,
                   const T& val0, const T& val1, const T& val2, const T& val3,
                   const T& pos);

/*! \brief Computes a specular direction. */
template <typename Vec3T>
Vec3T reflect(const Vec3T& dir, const Vec3T& normalDir);

/*! \brief Converts a value from radians to degrees. */
template <typename T>
T toDegree(T radian);

/*! \brief Converts a value from degrees to radians. */
template <typename T>
T toRadian(T degree);

/*! \brief Converts a coordinate system. */
template <typename SrcCoordSysT, typename DestCoordSysT>
void convertCoordinateSystem(float  srcAngle0,
                             float  srcAngle1,
                             float  srcAngle2,
                             float  srcAngle3,
                             float* destAngle0,
                             float* destAngle1,
                             float* destAngle2,
                             float* destAngle3);

/*! \brief Returns true if two sample sets have the same color model and wavelengths. */
template <typename T>
bool hasSameColor(const T& ss0, const T& ss1);

/*! \brief Converts from CIE-XYZ to sRGB. */
Vec3f xyzToSrgb(const Vec3f& xyz);

/*! \brief Converts from sRGB to CIE-XYZ. */
Vec3f srgbToXyz(const Vec3f& rgb);

/*! \brief Fixes a direction if the Z-component is negative. */
template <typename Vec3T>
void fixDownwardDir(Vec3T* dir);

/*! \brief Fits the angle with about the same value. */
template <typename T>
T fitAngle(T angle, T value);

/*
 * Implementation
 */

template <typename T>
inline T clamp(T value, T minValue, T maxValue)
{
    using std::min;
    using std::max;
    return max(minValue, min(maxValue, value));
}

template <typename T>
inline bool isEqual(T lhs, T rhs)
{
    using std::abs;
    using std::max;

    T tolerance = std::numeric_limits<T>::epsilon()
                * max(max(abs(lhs), abs(rhs)), static_cast<T>(1.0))
                * static_cast<T>(2.0);
    return (abs(lhs - rhs) <= tolerance);
}

template <typename T>
inline T lerp(const T& v0, const T& v1, float t)
{
    return v0 + (v1 - v0) * t;
}

template <typename T>
inline T smoothstep(const T& v0, const T& v1, const T& t)
{
    T coeff = clamp((t - v0) / (v1 - v0), 0.0f, 1.0f);
    return coeff * coeff * (3.0f - 2.0f * coeff);
}

template <typename T>
inline T smootherstep(const T& v0, const T& v1, const T& t)
{
    T coeff = clamp((t - v0) / (v1 - v0), 0.0f, 1.0f);
    return coeff * coeff * coeff * (coeff * (coeff * 6.0f - 15.0f) + 10.0f);
}

template <typename T>
inline T hermiteInterpolation3(const T& v0, const T& v1, float t)
{
    float coeff = smoothstep(0.0f, 1.0f, t);
    return lerp(v0, v1, coeff);
}

template <typename T>
inline T hermiteInterpolation5(const T& v0, const T& v1, float t)
{
    float coeff = smootherstep(0.0f, 1.0f, t);
    return lerp(v0, v1, coeff);
}

template <typename T>
inline T catmullRomSpline(const T& v0, const T& v1, const T& v2, const T& v3, float t)
{
    float t2 = t * t;
    float t3 = t2 * t;

    return static_cast<T>(((2.0 * v1) +
                           (-v0 + v2) * t +
                           (2.0 * v0 - 5.0 * v1 + 4.0 * v2 - v3) * t2 +
                           (-v0 + 3.0 * v1 - 3.0 * v2 + v3) * t3) * 0.5);
}

template <typename T>
inline T catmullRomSpline(const T& pos0, const T& pos1, const T& pos2, const T& pos3,
                          const T& val0, const T& val1, const T& val2, const T& val3,
                          const T& pos)
{
    Vec2 v0(pos0, val0);
    Vec2 v1(pos1, val1);
    Vec2 v2(pos2, val2);
    Vec2 v3(pos3, val3);
    CentripetalCatmullRomSpline ccrs(v0, v1, v2, v3);

    return static_cast<T>(ccrs.interpolateY(pos));
}

template <typename T>
inline bool hasSameColor(const T& ss0, const T& ss1)
{
    if (ss0.getColorModel() != ss1.getColorModel()) {
        std::cout
            << "[lb::hasSameColor] Color models do not match: "
            << ss0.getColorModel() << ", " << ss1.getColorModel()
            << std::endl;
        return false;
    }

    if (ss0.getNumWavelengths() != ss1.getNumWavelengths() ||
        !ss0.getWavelengths().isApprox(ss1.getWavelengths())) {
        std::cout
            << "[lb::hasSameColor] Wavelengths do not match: "
            << ss0.getWavelengths() << ", " << ss1.getWavelengths()
            << std::endl;
        return false;
    }

    return true;
}

template <typename Vec3T>
inline Vec3T reflect(const Vec3T& dir, const Vec3T& normalDir)
{
    typedef typename Vec3T::Scalar ScalarType;
    return static_cast<ScalarType>(2.0) * normalDir.dot(dir) * normalDir - dir;
}

template <typename T>
inline T toDegree(T radian)
{
    return radian / static_cast<T>(PI_F) * static_cast<T>(180.0);
}

template <typename T>
inline T toRadian(T degree)
{
    return degree / static_cast<T>(180.0) * static_cast<T>(PI_F);
}

template <typename SrcCoordSysT, typename DestCoordSysT>
inline void convertCoordinateSystem(float   srcAngle0,
                                    float   srcAngle1,
                                    float   srcAngle2,
                                    float   srcAngle3,
                                    float*  destAngle0,
                                    float*  destAngle1,
                                    float*  destAngle2,
                                    float*  destAngle3)
{
    Vec3 inDir, outDir;
    SrcCoordSysT::toXyz(srcAngle0, srcAngle1, srcAngle2, srcAngle3,
                        &inDir, &outDir);

    DestCoordSysT::fromXyz(inDir, outDir,
                           destAngle0, destAngle1, destAngle2, destAngle3);
}

inline Vec3f xyzToSrgb(const Vec3f& xyz)
{
    Eigen::Matrix3f mat;
    mat << CieData::XYZ_sRGB[0], CieData::XYZ_sRGB[1], CieData::XYZ_sRGB[2],
           CieData::XYZ_sRGB[3], CieData::XYZ_sRGB[4], CieData::XYZ_sRGB[5],
           CieData::XYZ_sRGB[6], CieData::XYZ_sRGB[7], CieData::XYZ_sRGB[8];
    return mat * xyz;
}

inline Vec3f srgbToXyz(const Vec3f& rgb)
{
    Eigen::Matrix3f mat;
    mat << CieData::sRGB_XYZ[0], CieData::sRGB_XYZ[1], CieData::sRGB_XYZ[2],
           CieData::sRGB_XYZ[3], CieData::sRGB_XYZ[4], CieData::sRGB_XYZ[5],
           CieData::sRGB_XYZ[6], CieData::sRGB_XYZ[7], CieData::sRGB_XYZ[8];
    return mat * rgb;
}

template <typename Vec3T>
inline void fixDownwardDir(Vec3T* dir)
{
    Vec3T& d = *dir;
    if (d[2] < 0.0) {
        d[2] = 0.0;
        if (d[0] == 0.0 && d[1] == 0.0) {
            d[0] = 1.0;
        }
        else {
            d.normalize();
        }
    }
}

template <typename T>
inline T fitAngle(T angle, T value)
{
    if (isEqual(angle, value)) {
        return value;
    }
    else {
        return angle;
    }
}

}  // namespace lb

#endif // LIBBSDF_UTILITY_H
