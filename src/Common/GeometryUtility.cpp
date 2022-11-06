// =================================================================== //
// Copyright (C) 2022 Kimura Ryo                                       //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#include <libbsdf/Common/GeometryUtility.h>

#include <libbsdf/Common/Global.h>

using namespace lb;

bool GeometryUtility::computeRayTriangleIntersection(const Vec3&   orig,
                                                     const Vec3&   dir,
                                                     const Vec3&   v0,
                                                     const Vec3&   v1,
                                                     const Vec3&   v2,
                                                     Vec3::Scalar* t,
                                                     Vec3::Scalar* u,
                                                     Vec3::Scalar* v)
{
    // Möller-Trumbore algorithm
    // https://www.scratchapixel.com/lessons/3d-basic-rendering/ray-tracing-rendering-a-triangle/moller-trumbore-ray-triangle-intersection

    using Scalar = Vec3::Scalar;

    Vec3   v0v1 = v1 - v0;
    Vec3   v0v2 = v2 - v0;
    Vec3   pVec = dir.cross(v0v2);
    Scalar det = v0v1.dot(pVec);

    // Ray and triangle are parallel.
    if (std::abs(det) < EPSILON_D) return false;

    Scalar invDet = 1 / det;

    Vec3 tVec = orig - v0;
    *u = tVec.dot(pVec) * invDet;
    if (*u < 0 || *u > 1) return false;

    Vec3 qVec = tVec.cross(v0v1);
    *v = dir.dot(qVec) * invDet;
    if (*v < 0 || *u + *v > 1) return false;

    *t = v0v2.dot(qVec) * invDet;

    return true;
}
