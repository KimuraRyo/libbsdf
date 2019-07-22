// =================================================================== //
// Copyright (C) 2014-2019 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

/*!
 * \file    Vector.h
 * \brief   The Vector.h header file includes the vector declarations and functions.
 */

#ifndef LIBBSDF_VECTOR_H
#define LIBBSDF_VECTOR_H

#include <Eigen/Core>
#include <Eigen/Geometry>

namespace lb {

typedef Eigen::Vector2f Vec2f;
typedef Eigen::Vector2d Vec2d;
typedef Eigen::Vector2i Vec2i;

typedef Eigen::Vector3f Vec3f;
typedef Eigen::Vector3d Vec3d;
typedef Eigen::Vector3i Vec3i;

typedef Eigen::Vector4f Vec4f;
typedef Eigen::Vector4d Vec4d;
typedef Eigen::Vector4i Vec4i;

#if defined(LIBBSDF_DOUBLE_PRECISION_VECTOR)
typedef Vec2d Vec2;
typedef Vec3d Vec3;
typedef Vec4d Vec4;
#else
typedef Vec2f Vec2;
typedef Vec3f Vec3;
typedef Vec4f Vec4;
#endif

/*! \brief Converts from a vector to lb::Vec3. */
template <typename Vec3T>
Vec3 toVec3(const Vec3T& vec3);

/*
 * Implementation
 */

template <typename Vec3T>
Vec3 toVec3(const Vec3T& vec3)
{
    return Vec3(static_cast<Vec3::Scalar>(vec3[0]),
                static_cast<Vec3::Scalar>(vec3[1]),
                static_cast<Vec3::Scalar>(vec3[2]));
}

} // namespace lb

#endif // LIBBSDF_VECTOR_H
