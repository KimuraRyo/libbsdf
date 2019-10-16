// =================================================================== //
// Copyright (C) 2019 Kimura Ryo                                       //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#include <libbsdf/Common/SolidAngle.h>

#include <libbsdf/Common/Log.h>

using namespace lb;

double SolidAngle::fromRectangleOnHemisphere(const Vec3&    v0,
                                             const Vec3&    v1,
                                             const Vec3&    v2,
                                             const Vec3&    v3,
                                             Vec3*          centroid)
{
    const Vec3::Scalar& z0 = v0.z();
    const Vec3::Scalar& z1 = v1.z();
    const Vec3::Scalar& z2 = v2.z();
    const Vec3::Scalar& z3 = v3.z();

    // All vertices are below the plane.
    if (z0 <= 0.0 && z1 <= 0.0 && z2 <= 0.0 && z3 <= 0.0) {
        return 0.0;
    }

    // All vertices are above the plane.
    else if (z0 >= 0.0 && z1 >= 0.0 && z2 >= 0.0 && z3 >= 0.0) {
        double area0 = fromTriangle(v0, v1, v2);
        double area1 = fromTriangle(v0, v2, v3);
        double area = area0 + area1;

        *centroid = (area0 * computeCentroid(v0, v1, v2) +
                     area1 * computeCentroid(v0, v2, v3)) / area;
        centroid->normalize();

        return area;
    }

    // One vertex is below the plane.
    else if (z0 < 0.0 && z1 >= 0.0 && z2 >= 0.0 && z3 >= 0.0) {
        Vec3 boundPos30 = computeBound(v3, v0);
        Vec3 boundPos10 = computeBound(v1, v0);

        double area0 = fromTriangle(v1, v2, v3);
        double area1 = fromTriangle(boundPos30, v1, v3);
        double area2 = fromTriangle(boundPos30, boundPos10, v1);
        double area = area0 + area1 + area2;

        *centroid = (area0 * computeCentroid(v1, v2, v3) +
                     area1 * computeCentroid(boundPos30, v1, v3) +
                     area2 * computeCentroid(boundPos30, boundPos10, v1)) / area;
        centroid->normalize();

        return area;
    }
    else if (z0 >= 0.0 && z1 < 0.0 && z2 >= 0.0 && z3 >= 0.0) {
        Vec3 boundPos01 = computeBound(v0, v1);
        Vec3 boundPos21 = computeBound(v2, v1);

        double area0 = fromTriangle(v2, v3, v0);
        double area1 = fromTriangle(boundPos01, v2, v0);
        double area2 = fromTriangle(boundPos01, boundPos21, v2);
        double area = area0 + area1 + area2;

        *centroid = (area0 * computeCentroid(v2, v3, v0) +
                     area1 * computeCentroid(boundPos01, v2, v0) +
                     area2 * computeCentroid(boundPos01, boundPos21, v2)) / area;
        centroid->normalize();

        return area;
    }
    else if (z0 >= 0.0 && z1 >= 0.0 && z2 < 0.0 && z3 >= 0.0) {
        Vec3 boundPos12 = computeBound(v1, v2);
        Vec3 boundPos32 = computeBound(v3, v2);

        double area0 = fromTriangle(v3, v0, v1);
        double area1 = fromTriangle(boundPos12, v3, v1);
        double area2 = fromTriangle(boundPos12, boundPos32, v3);
        double area = area0 + area1 + area2;

        *centroid = (area0 * computeCentroid(v3, v0, v1) +
                     area1 * computeCentroid(boundPos12, v3, v1) +
                     area2 * computeCentroid(boundPos12, boundPos32, v3)) / area;
        centroid->normalize();

        return area;
    }
    else if (z0 >= 0.0 && z1 >= 0.0 && z2 >= 0.0 && z3 < 0.0) {
        Vec3 boundPos23 = computeBound(v2, v3);
        Vec3 boundPos03 = computeBound(v0, v3);

        double area0 = fromTriangle(v0, v1, v2);
        double area1 = fromTriangle(boundPos23, v0, v2);
        double area2 = fromTriangle(boundPos23, boundPos03, v0);
        double area = area0 + area1 + area2;

        *centroid = (area0 * computeCentroid(v0, v1, v2) +
                     area1 * computeCentroid(boundPos23, v0, v2) +
                     area2 * computeCentroid(boundPos23, boundPos03, v0)) / area;
        centroid->normalize();

        return area;
    }

    // Two vertices are below the plane.
    else if (z0 > 0.0 && z1 > 0.0 && z2 < 0.0 && z3 < 0.0) {
        Vec3 boundPos03 = computeBound(v0, v3);
        Vec3 boundPos12 = computeBound(v1, v2);

        double area0 = fromTriangle(v0, v1, boundPos12);
        double area1 = fromTriangle(v0, boundPos12, boundPos03);
        double area = area0 + area1;

        *centroid = (area0 * computeCentroid(v0, v1, boundPos12) +
                     area1 * computeCentroid(v0, boundPos12, boundPos03)) / area;
        centroid->normalize();

        return area;
    }
    else if (z0 < 0.0 && z1 > 0.0 && z2 > 0.0 && z3 < 0.0) {
        Vec3 boundPos10 = computeBound(v1, v0);
        Vec3 boundPos23 = computeBound(v2, v3);

        double area0 = fromTriangle(v1, v2, boundPos23);
        double area1 = fromTriangle(v1, boundPos23, boundPos10);
        double area = area0 + area1;

        *centroid = (area0 * computeCentroid(v1, v2, boundPos23) +
                     area1 * computeCentroid(v1, boundPos23, boundPos10)) / area;
        centroid->normalize();

        return area;
    }
    else if (z0 < 0.0 && z1 < 0.0 && z2 > 0.0 && z3 > 0.0) {
        Vec3 boundPos21 = computeBound(v2, v1);
        Vec3 boundPos30 = computeBound(v3, v0);

        double area0 = fromTriangle(v2, v3, boundPos30);
        double area1 = fromTriangle(v2, boundPos30, boundPos21);
        double area = area0 + area1;

        *centroid = (area0 * computeCentroid(v2, v3, boundPos30) +
                     area1 * computeCentroid(v2, boundPos30, boundPos21)) / area;
        centroid->normalize();

        return area;
    }
    else if (z0 > 0.0 && z1 < 0.0 && z2 < 0.0 && z3 > 0.0) {
        Vec3 boundPos32 = computeBound(v3, v2);
        Vec3 boundPos01 = computeBound(v0, v1);

        double area0 = fromTriangle(v3, v0, boundPos01);
        double area1 = fromTriangle(v3, boundPos01, boundPos32);
        double area = area0 + area1;

        *centroid = (area0 * computeCentroid(v3, v0, boundPos01) +
                     area1 * computeCentroid(v3, boundPos01, boundPos32)) / area;
        centroid->normalize();

        return area;
    }

    // Three vertices are below the plane.
    else if (z0 > 0.0 && z1 <= 0.0 && z2 <= 0.0 && z3 <= 0.0) {
        Vec3 boundPos01 = computeBound(v0, v1);
        Vec3 boundPos03 = computeBound(v0, v3);

        *centroid = computeCentroid(v0, boundPos01, boundPos03);
        centroid->normalize();

        return fromTriangle(v0, boundPos01, boundPos03);
    }
    else if (z0 <= 0.0 && z1 > 0.0 && z2 <= 0.0 && z3 <= 0.0) {
        Vec3 boundPos12 = computeBound(v1, v2);
        Vec3 boundPos10 = computeBound(v1, v0);

        *centroid = computeCentroid(v1, boundPos12, boundPos10);
        centroid->normalize();

        return fromTriangle(v1, boundPos12, boundPos10);
    }
    else if (z0 <= 0.0 && z1 <= 0.0 && z2 > 0.0 && z3 <= 0.0) {
        Vec3 boundPos23 = computeBound(v2, v3);
        Vec3 boundPos21 = computeBound(v2, v1);

        *centroid = computeCentroid(v2, boundPos23, boundPos21);
        centroid->normalize();

        return fromTriangle(v2, boundPos23, boundPos21);
    }
    else if (z0 <= 0.0 && z1 <= 0.0 && z2 <= 0.0 && z3 > 0.0) {
        Vec3 boundPos30 = computeBound(v3, v0);
        Vec3 boundPos32 = computeBound(v3, v2);

        *centroid = computeCentroid(v3, boundPos30, boundPos32);
        centroid->normalize();

        return fromTriangle(v3, boundPos30, boundPos32);
    }

    else {
        lbError << "[SolidAngle::fromRectangleOnHemisphere] Failed to compute a solid angle.";
        return 0.0;
    }
}
