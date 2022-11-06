// =================================================================== //
// Copyright (C) 2022 Kimura Ryo                                       //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#ifndef LIBBSDF_GEOMETRY_UTILITY_H
#define LIBBSDF_GEOMETRY_UTILITY_H

#include <libbsdf/Common/Vector.h>

namespace lb {

/*!
 * \class   GeometryUtility
 * \brief   The GeometryUtility class provides utility functions for geometry.
 */
class GeometryUtility
{
public:
    /*!
     * Computes the ray-triangle intersection.
     * 
     * \note    The third value of the barycentric coordinates can be calculated from \a u and \a v.
     *          w = 1 - u - v
     * 
     * \param t The distance from the ray's origin to the intersection.
     * \param u The barycentric coordinate.
     * \param v The barycentric coordinate.
     */
    static bool computeRayTriangleIntersection(const Vec3&   orig,
                                               const Vec3&   dir,
                                               const Vec3&   v0,
                                               const Vec3&   v1,
                                               const Vec3&   v2,
                                               Vec3::Scalar* t,
                                               Vec3::Scalar* u,
                                               Vec3::Scalar* v);

    /*! Computes the barycentric coordinate of \a p in a triangle abc. */
    template <typename Vec2T, typename Vec3T>
    static Vec3T computeBarycentricCoord(const Vec2T& p, const Vec2T& a, const Vec2T& b, const Vec2T& c);
};

template <typename Vec2T, typename Vec3T>
Vec3T GeometryUtility::computeBarycentricCoord(const Vec2T& p,
                                               const Vec2T& a,
                                               const Vec2T& b,
                                               const Vec2T& c)
{
    using V2Scalar = typename Vec2T::Scalar;
    using V3Scalar = typename Vec3T::Scalar;

    Vec2T v0 = b - a;
    Vec2T v1 = c - a;
    Vec2T v2 = p - a;

    V2Scalar denom = v0.x() * v1.y() - v1.x() * v0.y();

    V3Scalar v = static_cast<V3Scalar>((v2.x() * v1.y() - v1.x() * v2.y()) / denom);
    V3Scalar w = static_cast<V3Scalar>((v0.x() * v2.y() - v2.x() * v0.y()) / denom);
    V3Scalar u = V3Scalar(1) - v - w;

    return Vec3T(u, v, w);
}

} // namespace lb

#endif // LIBBSDF_GEOMETRY_UTILITY_H
