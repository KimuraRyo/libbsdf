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

using Vec2f = Eigen::Vector2f;
using Vec2d = Eigen::Vector2d;
using Vec2i = Eigen::Vector2i;

using Vec3f = Eigen::Vector3f;
using Vec3d = Eigen::Vector3d;
using Vec3i = Eigen::Vector3i;

using Vec4f = Eigen::Vector4f;
using Vec4d = Eigen::Vector4d;
using Vec4i = Eigen::Vector4i;

#if defined(LIBBSDF_DOUBLE_PRECISION_VECTOR)
using Vec2 = Vec2d;
using Vec3 = Vec3d;
using Vec4 = Vec4d;
#else
using Vec2 = Vec2f;
using Vec3 = Vec3f;
using Vec4 = Vec4f;
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
