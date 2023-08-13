// =================================================================== //
// Copyright (C) 2020-2023 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#ifndef LIBBSDF_MONOTONE_CUBIC_INTERPOLATION_H
#define LIBBSDF_MONOTONE_CUBIC_INTERPOLATION_H

#include <cassert>

#include <libbsdf/Common/Vector.h>

namespace lb {

/*!
 * \class   MonotoneCubicInterpolation
 * \brief   The MonotoneCubicInterpolation class provides the functions for monotone cubic interpolation.
 */
class MonotoneCubicInterpolation
{
public:
    /*! Computes an interpolated value using monotone cubic interpolation at \a x in [\a v1.x(),\a v2.x()]. */
    template <typename Vec2T>
    static typename Vec2T::Scalar interpolate(const Vec2T&                  v0,
                                              const Vec2T&                  v1,
                                              const Vec2T&                  v2,
                                              const Vec2T&                  v3,
                                              const typename Vec2T::Scalar& x);

    /*!
     * Computes an interpolated array using monotone cubic interpolation at \a pos in [\a pos1,\a pos2].
     *
     * \return Interpolated array.
     */
    template <typename ArrayT>
    static ArrayT interpolate(double        pos0,
                              double        pos1,
                              double        pos2,
                              double        pos3,
                              const ArrayT& arr0,
                              const ArrayT& arr1,
                              const ArrayT& arr2,
                              const ArrayT& arr3,
                              double        pos);

private:
    /*!
     * Computes the slope at \a currV.
     * See F. N. Fritsch and R. E. Carlson. 1980. "Monotone Piecewise Cubic Interpolation".
     */
    template <typename Vec2T>
    static typename Vec2T::Scalar computeSlope(const Vec2T& prevV,
                                               const Vec2T& currV,
                                               const Vec2T& nextV);

    /*!
     * Computes the slope at \a currArr.
     * See F. N. Fritsch and R. E. Carlson. 1980. "Monotone Piecewise Cubic Interpolation".
     */
    template <typename ArrayT>
    static ArrayT computeSlope(double        prevPos,
                               double        currPos,
                               double        nextPos,
                               const ArrayT& prevArr,
                               const ArrayT& currArr,
                               const ArrayT& nextArr);
};

template <typename Vec2T>
typename Vec2T::Scalar MonotoneCubicInterpolation::interpolate(const Vec2T&                  v0,
                                                               const Vec2T&                  v1,
                                                               const Vec2T&                  v2,
                                                               const Vec2T&                  v3,
                                                               const typename Vec2T::Scalar& x)
{
    // Reference:
    // https://math.stackexchange.com/questions/4082/equation-of-a-curve-given-3-points-and-additional-constant-requirements#4104
    // https://math.stackexchange.com/questions/45218/implementation-of-monotone-cubic-interpolation#51412

    using Scalar = typename Vec2T::Scalar;

    if (v1.x() == v2.x()) {
        return (v1.y() + v2.y()) * Scalar(0.5);
    }

    Scalar slope1 = computeSlope(v0, v1, v2);
    Scalar slope2 = computeSlope(v1, v2, v3);

    Vec2T diff = v2 - v1;
    Scalar distX = diff.x();
    Scalar slope = diff.y() / distX;
    Scalar posX = x - v1.x();

    return v1.y() +
           slope1 * posX +
           (3 * slope - 2 * slope1 - slope2) / distX * posX * posX +
           (slope1 + slope2 - 2 * slope) / (distX * distX) * posX * posX * posX;
}

template <typename ArrayT>
ArrayT MonotoneCubicInterpolation::interpolate(double        pos0,
                                               double        pos1,
                                               double        pos2,
                                               double        pos3,
                                               const ArrayT& arr0,
                                               const ArrayT& arr1,
                                               const ArrayT& arr2,
                                               const ArrayT& arr3,
                                               double        pos)
{
    assert(arr0.size() == arr1.size() &&
           arr1.size() == arr2.size() &&
           arr2.size() == arr3.size());

    using Scalar = typename ArrayT::Scalar;

    if (pos1 == pos2) {
        return (arr1 + arr2) * Scalar(0.5);
    }

    ArrayT slope1 = computeSlope(pos0, pos1, pos2, arr0, arr1, arr2);
    ArrayT slope2 = computeSlope(pos1, pos2, pos3, arr1, arr2, arr3);

    Scalar distX = static_cast<Scalar>(pos2 - pos1);
    ArrayT slope = (arr2 - arr1) / distX;
    Scalar posX = static_cast<Scalar>(pos - pos1);

    ArrayT arr = arr1 +
                 slope1 * posX +
                 (3 * slope - 2 * slope1 - slope2) / distX * posX * posX +
                 (slope1 + slope2 - 2 * slope) / (distX * distX) * posX * posX * posX;

    assert(arr.allFinite());

    return arr;
}

template <typename Vec2T>
typename Vec2T::Scalar MonotoneCubicInterpolation::computeSlope(const Vec2T& prevV,
                                                                const Vec2T& currV,
                                                                const Vec2T& nextV)
{
    using Scalar = typename Vec2T::Scalar;

    Vec2T prevDiff = currV - prevV;
    Vec2T nextDiff = nextV - currV;
    Scalar prevDistX = prevDiff.x();
    Scalar nextDistX = nextDiff.x();

    Scalar prevSlope = (prevDistX == Scalar(0)) ? Scalar(0) : prevDiff.y() / prevDistX;
    Scalar nextSlope = (nextDistX == Scalar(0)) ? Scalar(0) : nextDiff.y() / nextDistX;

    if (sign(prevSlope) == sign(nextSlope)) {
        return 3 * (prevDistX + nextDistX) / ((2 * nextDistX + prevDistX) / prevSlope +
                                              (nextDistX + 2 * prevDistX) / nextSlope);
    }
    else {
        return Scalar(0);
    }
}

template <typename ArrayT>
ArrayT MonotoneCubicInterpolation::computeSlope(double        prevPos,
                                                double        currPos,
                                                double        nextPos,
                                                const ArrayT& prevArr,
                                                const ArrayT& currArr,
                                                const ArrayT& nextArr)
{
    using Scalar = typename ArrayT::Scalar;

    ArrayT prevDiff = currArr - prevArr;
    ArrayT nextDiff = nextArr - currArr;
    Scalar prevDistX = static_cast<Scalar>(currPos - prevPos);
    Scalar nextDistX = static_cast<Scalar>(nextPos - currPos);

    ArrayT prevSlope;
    if (prevDistX == Scalar(0)) {
        prevSlope = ArrayT::Zero(prevDiff.size());
    }
    else {
        prevSlope = prevDiff / prevDistX;
    }

    ArrayT nextSlope;
    if (nextDistX == Scalar(0)) {
        nextSlope = ArrayT::Zero(nextDiff.size());
    }
    else {
        nextSlope = nextDiff / nextDistX;
    }

    ArrayT slope(currArr.size());

    for (int i = 0; i < currArr.size(); ++i) {
        Scalar ps = prevSlope[i];
        Scalar ns = nextSlope[i];

        if (ps != Scalar(0) &&
            ns != Scalar(0) &&
            sign(ps) == sign(ns)) {
            slope[i] = 3 * (prevDistX + nextDistX) / ((2 * nextDistX + prevDistX) / ps +
                                                      (nextDistX + 2 * prevDistX) / ns);
        }
        else {
            slope[i] = Scalar(0);
        }
    }

    return slope;
}

} // namespace lb

#endif // LIBBSDF_MONOTONE_CUBIC_INTERPOLATION_H
