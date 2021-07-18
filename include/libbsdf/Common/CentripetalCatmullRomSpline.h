// =================================================================== //
// Copyright (C) 2015-2020 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#ifndef LIBBSDF_CENTRIPETAL_CATMULL_ROM_SPLINE_H
#define LIBBSDF_CENTRIPETAL_CATMULL_ROM_SPLINE_H

#include <cassert>

#include <libbsdf/Common/Vector.h>

namespace lb {

/*!
 * \class   CentripetalCatmullRomSpline
 * \brief   The CentripetalCatmullRomSpline class provides the interpolation using centripetal Catmull-Rom spline.
 */
class CentripetalCatmullRomSpline
{
public:
    CentripetalCatmullRomSpline();

    CentripetalCatmullRomSpline(const Vec2& pos0,
                                const Vec2& pos1,
                                const Vec2& pos2,
                                const Vec2& pos3);

    /*! Initializes parameters. */
    void initialize(const Vec2& pos0,
                    const Vec2& pos1,
                    const Vec2& pos2,
                    const Vec2& pos3);

    /*! Evaluates the spline at \a t in [0,1]. */
    Vec2 evaluate(const Vec2::Scalar& t);

    /*! Interpolates the Y-component of position using X-component. */
    Vec2::Scalar interpolateY(const Vec2::Scalar& x);

    /*! Computes an interpolated value using centripetal Catmull-Rom spline at \a x in [\a v1.x(),\a v2.x()]. */
    template <typename Vec2T>
    static typename Vec2T::Scalar interpolate(const Vec2T&                  v0,
                                              const Vec2T&                  v1,
                                              const Vec2T&                  v2,
                                              const Vec2T&                  v3,
                                              const typename Vec2T::Scalar& x);

    /*!
     * Computes an interpolated array using centripetal Catmull-Rom spline at \a pos in [\a pos1,\a pos2].
     *
     * \return Interpolated array.
     */
    template <typename ArrayT>
    static ArrayT interpolate(const typename ArrayT::Scalar&    pos0,
                              const typename ArrayT::Scalar&    pos1,
                              const typename ArrayT::Scalar&    pos2,
                              const typename ArrayT::Scalar&    pos3,
                              const ArrayT&                     arr0,
                              const ArrayT&                     arr1,
                              const ArrayT&                     arr2,
                              const ArrayT&                     arr3,
                              const typename ArrayT::Scalar&    pos);

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

template <typename Vec2T>
typename Vec2T::Scalar CentripetalCatmullRomSpline::interpolate(const Vec2T&                    v0,
                                                                const Vec2T&                    v1,
                                                                const Vec2T&                    v2,
                                                                const Vec2T&                    v3,
                                                                const typename Vec2T::Scalar&   x)
{
    CentripetalCatmullRomSpline ccrs(v0.template cast<Vec2::Scalar>(),
                                     v1.template cast<Vec2::Scalar>(),
                                     v2.template cast<Vec2::Scalar>(),
                                     v3.template cast<Vec2::Scalar>());

    using Scalar = typename Vec2T::Scalar;
    return static_cast<Scalar>(ccrs.interpolateY(static_cast<double>(x)));
}

template <typename ArrayT>
ArrayT CentripetalCatmullRomSpline::interpolate(const typename ArrayT::Scalar&  pos0,
                                                const typename ArrayT::Scalar&  pos1,
                                                const typename ArrayT::Scalar&  pos2,
                                                const typename ArrayT::Scalar&  pos3,
                                                const ArrayT&                   arr0,
                                                const ArrayT&                   arr1,
                                                const ArrayT&                   arr2,
                                                const ArrayT&                   arr3,
                                                const typename ArrayT::Scalar&  pos)
{
    assert(arr0.size() == arr1.size() &&
           arr1.size() == arr2.size() &&
           arr2.size() == arr3.size());

    ArrayT arr;
    arr.resize(arr0.size());

    CentripetalCatmullRomSpline ccrs;
    #pragma omp parallel for private(ccrs)
    for (int i = 0; i < arr.size(); ++i) {
        ccrs.initialize(Vec2(pos0, arr0[i]),
                        Vec2(pos1, arr1[i]),
                        Vec2(pos2, arr2[i]),
                        Vec2(pos3, arr3[i]));

        using Scalar = typename ArrayT::Scalar;
        arr[i] = static_cast<Scalar>(ccrs.interpolateY(pos));
    }

    return arr;
}

} // namespace lb

#endif // LIBBSDF_CENTRIPETAL_CATMULL_ROM_SPLINE_H
