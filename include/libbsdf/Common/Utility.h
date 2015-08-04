// =================================================================== //
// Copyright (C) 2014-2015 Kimura Ryo                                  //
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

#include <libbsdf/Common/Global.h>

namespace lb {

/*! \brief Clamps a value between a minimum and maximum value. */
template <typename T>
T clamp(T value, T minValue, T maxValue);

/*! \brief Returns true if two values are nearly equal. */
template <typename T>
bool isEqual(T lhs, T rhs);

/*! \brief Computes linearly-interpolated values. */
template <typename T>
T lerp(const T& lhs, const T& rhs, float weight);

/*! \brief Computes a specular direction. */
template <typename Vec3T>
Vec3T reflect(const Vec3T& inDir, const Vec3T& normalDir);

/*! \brief Converts a value from radian to degree. */
template <typename T>
T toDegree(T radian);

/*! \brief Converts a value from degree to radian. */
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

/*! \brief Fixes a direction if the Z-component is negative. */
template <typename Vec3T>
void fixDownwardDir(Vec3T& dir);

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
    return (abs(lhs - rhs) <= std::numeric_limits<T>::epsilon() * abs(lhs + rhs));
}

template <typename T>
inline T lerp(const T& lhs, const T& rhs, float weight)
{
    return lhs + weight * (rhs - lhs);
}

template <typename Vec3T>
inline Vec3T reflect(const Vec3T& inDir, const Vec3T& normalDir)
{
    typedef typename Vec3T::Scalar ScalarType;
    return static_cast<ScalarType>(2.0) * normalDir.dot(inDir) * normalDir - inDir;
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

template <typename Vec3T>
inline void fixDownwardDir(Vec3T& dir)
{
    if (dir[2] < 0.0) {
        dir[2] = 0.0;
        if (dir[0] == 0.0 && dir[1] == 0.0) {
            dir[0] = 1.0;
        }
        else {
            dir.normalize();
        }
    }
}

}  // namespace lb

#endif // LIBBSDF_UTILITY_H
