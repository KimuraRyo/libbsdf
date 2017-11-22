// =================================================================== //
// Copyright (C) 2015-2017 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#include <libbsdf/Common/CentripetalCatmullRomSpline.h>

#include <cassert>

#include <libbsdf/Common/Utility.h>

using namespace lb;

CentripetalCatmullRomSpline::CentripetalCatmullRomSpline() {}

CentripetalCatmullRomSpline::CentripetalCatmullRomSpline(const Vec2& pos0,
                                                         const Vec2& pos1,
                                                         const Vec2& pos2,
                                                         const Vec2& pos3)
{
    initialize(pos0, pos1, pos2, pos3);
}

void CentripetalCatmullRomSpline::initialize(const Vec2& pos0,
                                             const Vec2& pos1,
                                             const Vec2& pos2,
                                             const Vec2& pos3)
{
    using std::max;
    using std::sqrt;

    pos0_ = pos0;
    pos1_ = pos1;
    pos2_ = pos2;
    pos3_ = pos3;

    // Compute disntaces.
    Vec2::Scalar dist0 = (pos1 - pos0).norm();
    Vec2::Scalar dist1 = (pos2 - pos1).norm();
    Vec2::Scalar dist2 = (pos3 - pos2).norm();

    // Modify distances for centripetal Catmull-Rom spline.
    dist0 = max(sqrt(dist0), std::numeric_limits<Vec2::Scalar>::epsilon());
    dist1 = max(sqrt(dist1), std::numeric_limits<Vec2::Scalar>::epsilon());
    dist2 = max(sqrt(dist2), std::numeric_limits<Vec2::Scalar>::epsilon());

    // Compute tangents.
    Vec2 tan1 = (pos1 - pos0) / dist0 - (pos2 - pos0) / (dist0 + dist1) + (pos2 - pos1) / dist1;
    Vec2 tan2 = (pos2 - pos1) / dist1 - (pos3 - pos1) / (dist1 + dist2) + (pos3 - pos2) / dist2;
    tan1 *= dist1;
    tan2 *= dist1;

    computeCoefficients(pos1, pos2, tan1, tan2);
}

Vec2::Scalar CentripetalCatmullRomSpline::interpolateY(const Vec2::Scalar& x)
{
    assert((x > pos1_.x() - EPSILON_F * 10.0f && x < pos2_.x() + EPSILON_F * 10.0f) ||
           (x < pos1_.x() + EPSILON_F * 10.0f && x > pos2_.x() - EPSILON_F * 10.0f));

    Vec2::Scalar t = (x - pos1_.x()) / (pos2_.x() - pos1_.x());

    using std::abs;
    if (abs(pos1_.x() - x) < EPSILON_F * 10.0f) return pos1_.y();
    if (abs(pos2_.x() - x) < EPSILON_F * 10.0f) return pos2_.y();

    Vec2 currentPos = evaluate(t);

    Vec2::Scalar minT = 0.0;
    Vec2::Scalar maxT = 1.0;

    if (pos1_.x() < pos2_.x()) {
        while (!isEqual(currentPos.x(), x) &&
               t != minT &&
               t != maxT) {
            if (currentPos.x() > x) {
                maxT = t;
                t = (minT + t) * 0.5;
            }
            else {
                minT = t;
                t = (t + maxT) * 0.5;
            }

            currentPos = evaluate(t);
        }
    }
    else {
        while (!isEqual(currentPos.x(), x) &&
               t != minT &&
               t != maxT) {
            if (currentPos.x() < x) {
                maxT = t;
                t = (minT + t) * 0.5;
            }
            else {
                minT = t;
                t = (t + maxT) * 0.5;
            }

            currentPos = evaluate(t);
        }
    }

    return currentPos.y();
}
