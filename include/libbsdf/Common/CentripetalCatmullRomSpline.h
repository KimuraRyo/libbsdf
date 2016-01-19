// =================================================================== //
// Copyright (C) 2015 Kimura Ryo                                       //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#ifndef LIBBSDF_CENTRIPETAL_CATMULL_ROM_SPLINE_H
#define LIBBSDF_CENTRIPETAL_CATMULL_ROM_SPLINE_H

#include <libbsdf/Common/Vector.h>

namespace lb {

/*!
 * \class   CentripetalCatmullRomSpline
 * \brief   The CentripetalCatmullRomSpline class provides the interpolation using
 *          centripetal Catmull-Rom spline.
 */
class CentripetalCatmullRomSpline
{
public:
    CentripetalCatmullRomSpline(const Vec2& pos0,
                                const Vec2& pos1,
                                const Vec2& pos2,
                                const Vec2& pos3);

    /*! Evaluates the spline at \a t in [0,1]. */
    Vec2 evaluate(const Vec2::Scalar& t);

    /*! Interpolates the Y-component of position using X-component. */
    Vec2::Scalar interpolateY(const Vec2::Scalar& x);

private:
    /*! Computes the coefficients of cubic Hermite spline using positions and tangents. */
    void computeCoefficients(const Vec2& pos1,
                             const Vec2& pos2,
                             const Vec2& tan1,
                             const Vec2& tan2);

    Vec2 pos0_, pos1_, pos2_, pos3_;
    Vec2 coeff0_, coeff1_, coeff2_, coeff3_;
};

inline Vec2 CentripetalCatmullRomSpline::evaluate(const Vec2::Scalar& t)
{
    Vec2::Scalar t2 = t * t;
    Vec2::Scalar t3 = t2 * t;
    return coeff0_ + coeff1_ * t + coeff2_ * t2 + coeff3_ * t3;
}

inline void CentripetalCatmullRomSpline::computeCoefficients(const Vec2& pos1,
                                                             const Vec2& pos2,
                                                             const Vec2& tan1,
                                                             const Vec2& tan2)
{
    coeff0_ = pos1;
    coeff1_ = tan1;
    coeff2_ = -3.0 * pos1 + 3.0 * pos2 - 2.0 * tan1 - tan2;
    coeff3_ = 2.0 * pos1 - 2.0 * pos2 + tan1 + tan2;
}

} // namespace lb

#endif // LIBBSDF_CENTRIPETAL_CATMULL_ROM_SPLINE_H
