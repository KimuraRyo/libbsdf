// =================================================================== //
// Copyright (C) 2022 Kimura Ryo                                       //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#ifndef LIBBSDF_STEREOGRAPHIC_PROJECTION_H
#define LIBBSDF_STEREOGRAPHIC_PROJECTION_H

#include <libbsdf/Common/Vector.h>

namespace lb {

/*!
 * \struct  StereographicProjection
 * \brief   The StereographicProjection struct provides the functions of the stereographic projection.
 *
 * The formulation is based on the projection from the south pole (z = -1) onto the plane through the equator (z = 0).
 * Points in the northern hemisphere (z = [0, 1]) are mapped inside the unit circle.
 * 
 */
struct StereographicProjection
{
    /*! Converts Cartesian coordinates on the unit sphere to points on the plane. */
    static Vec2 fromXyz(Vec3 xyz);

    /*! Converts points on the plane to Cartesian coordinates on the unit sphere. */
    static Vec3 toXyz(Vec2 xy);
};

inline Vec2 StereographicProjection::fromXyz(Vec3 xyz)
{
    return Vec2(xyz.x(), xyz.y()) / (1 + xyz.z());
}

inline Vec3 StereographicProjection::toXyz(Vec2 xy)
{
    Vec3::Scalar x2 = xy.x() * xy.x();
    Vec3::Scalar y2 = xy.y() * xy.y();
    return Vec3(2 * xy.x(), 2 * xy.y(), 1 - x2 - y2) / (1 + x2 + y2);
}

} // namespace lb

#endif // LIBBSDF_STEREOGRAPHIC_PROJECTION_H
