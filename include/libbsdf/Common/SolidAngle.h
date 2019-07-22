// =================================================================== //
// Copyright (C) 2019 Kimura Ryo                                       //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#ifndef SOLID_ANGLE_UTILITY_H
#define SOLID_ANGLE_UTILITY_H

#include <cmath>

#include <libbsdf/Common/Vector.h>

namespace lb {

/*!
 * \class   SolidAngle
 * \brief   The SolidAngle class provides utility functions for solid angle.
 */
class SolidAngle
{
public:
    /*! Computes the solid angle of a rectangle specified by polar and azimuthal angles. */
    template <typename T>
    static T fromRectangle(const T& theta0, const T& theta1, const T& phi0, const T& phi1);

    /*! 
     * Computes the solid angle of a rectangle specified by four unit vectors.
     * If the part of a rectangle is below the plane (z=0), the upper region is validated.
     * Vertices of a rectangle:
     *   v0--v1
     *   |   |
     *   v3--v2
     * 
     * \param   The assigned centroid of a rectangle.
     *          If the part of a rectangle is below the plane, the upper part is taken into account.
     *          If all vectors are below the plane, result is not defined.
     * \return  The solid angle above the plane.
     */
    static double fromRectangleOnHemisphere(const Vec3& v0,
                                            const Vec3& v1,
                                            const Vec3& v2,
                                            const Vec3& v3,
                                            Vec3*       centroid);

    /*! Computes the solid angle of a triangle specified by three unit vectors. */
    static double fromTriangle(const Vec3& v0, const Vec3& v1, const Vec3& v2);

private:
    /*!
     * Computes the centroid of a planar triangle.
     * This is different from the centroid of a spherial triangle, but similar for a small triangle.
     */
    static Vec3 computeCentroid(const Vec3& v0, const Vec3& v1, const Vec3& v2);

    /*! Computes the intersection of the plane (z=0) and a line segment with \a v0 and \a v1. */
    static Vec3 computeBound(const Vec3& v0, const Vec3& v1);
};

template <typename T>
inline T SolidAngle::fromRectangle(const T& theta0, const T& theta1, const T& phi0, const T& phi1)
{
    using std::abs;
    using std::cos;

    return abs((cos(theta0) - cos(theta1)) * (phi0 - phi1));
}

inline double SolidAngle::fromTriangle(const Vec3& v0, const Vec3& v1, const Vec3& v2)
{
    using std::abs;
    using std::atan;

    return abs(atan(v0.dot(v1.cross(v2)) / (1.0 + v0.dot(v1) + v1.dot(v2) + v2.dot(v0)))) * 2.0;
}

inline Vec3 SolidAngle::computeCentroid(const Vec3& v0, const Vec3& v1, const Vec3& v2)
{
    return ((v0 + v1 + v2) / 3.0).normalized();
}

inline Vec3 SolidAngle::computeBound(const Vec3& v0, const Vec3& v1)
{
    return (v1 + (-v1[2] / (v0[2] - v1[2])) * (v0 - v1)).normalized();
}

} // namespace lb

#endif // SOLID_ANGLE_UTILITY_H
