// =================================================================== //
// Copyright (C) 2022 Kimura Ryo                                       //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#include <libbsdf/Brdf/ScatteredSampleSet2D.h>

#include <libbsdf/Common/GeometryUtility.h>
#include <libbsdf/Common/Utility.h>

using namespace lb;

ScatteredSampleSet2D::~ScatteredSampleSet2D() {}

bool ScatteredSampleSet2D::constructData(bool extrapolationOfSamplesOnUnitSphere)
{
    if (sampleMap_.size() < 3) {
        lbWarn << "[ScatteredSampleSet2D::constructData] The number of samples aren't enough. Size: "
               << sampleMap_.size();
        return false;
    }

    dt_.initialize(sampleMap_.size());
    spectra_.reserve(sampleMap_.size());

    for (const auto& s : sampleMap_) {
        dt_.addVertex(s.first);
        spectra_.push_back(s.second);
    }

    dt_.computeDelaunayTriangles();

    addExtrapolatedSamplesOnUnitSphere();

    return true;
}

Spectrum ScatteredSampleSet2D::getSpectrum(const Vec3& dir) const
{
    for (size_t i = 0; i < dt_.getNumTriangles(); ++i) {
        Vec3i trg(dt_.getTriangle(i));

        size_t index0 = trg[0];
        size_t index1 = trg[1];
        size_t index2 = trg[2];

        Vec3 vert0(StereographicProjection::toXyz(dt_.getVertex(index0)));
        Vec3 vert1(StereographicProjection::toXyz(dt_.getVertex(index1)));
        Vec3 vert2(StereographicProjection::toXyz(dt_.getVertex(index2)));

        // Make the triangle slightly larger to detect the intersection on edge or vertex.
        constexpr Vec3::Scalar eps = EPSILON * 10;
        Vec3 center = (vert0 + vert1 + vert2) / 3;
        vert0 += (vert0 - center) * eps;
        vert1 += (vert1 - center) * eps;
        vert2 += (vert2 - center) * eps;

        Spectrum sp0(spectra_.at(index0));
        Spectrum sp1(spectra_.at(index1));
        Spectrum sp2(spectra_.at(index2));

        Vec3         orig(0, 0, 0);
        Vec3::Scalar t, u, v;
        if (GeometryUtility::computeRayTriangleIntersection(orig, dir, vert0, vert1, vert2, &t, &u,
                                                            &v)) {
            return (1 - u - v) * sp0 + u * sp1 + v * sp2;
        }
    }

    return lb::Spectrum();
}

Vec2 ScatteredSampleSet2D::computeIntersectionOnUnitCircle(const Vec2& pos, const Vec2& dir)
{
    Vec2         u = -pos;
    Vec2         u1 = u.dot(dir) * dir;
    Vec2         u2 = u - u1;
    Vec2::Scalar d = u2.norm();
    Vec2::Scalar m = std::sqrt(1 - d * d);

    Vec2 p1 = pos + u1 + m * dir;
    Vec2 p2 = pos + u1 - m * dir;

    if (dir.dot((p1 - pos).normalized()) > 0) {
        return p1;
    }
    else {
        return p2;
    }
}

void ScatteredSampleSet2D::addExtrapolatedSamplesOnUnitSphere()
{
    const std::vector<Vec2>& convexHull = dt_.getConvexHull();
    for (size_t i = 0; i < convexHull.size(); ++i) {
        size_t backI = convexHull.size() - 1;
        size_t prevI = (i != 0) ? i - 1 : backI;
        size_t nextI = (i != backI) ? i + 1 : 0;

        Vec2 pos = convexHull.at(i);

        Vec2 prevEdge = pos - convexHull.at(prevI);
        Vec2 nextEdge = pos - convexHull.at(nextI);
        Vec2 dir = (prevEdge.normalized() + nextEdge.normalized()).normalized();

        Vec2 v = computeIntersectionOnUnitCircle(pos, dir);

        constexpr Vec2::Scalar espDist = EPSILON_F;

        Vec2 insideConvexHullPos = pos - espDist * dir;

        Vec3 onConvexHullDir = StereographicProjection::toXyz(pos);
        Vec3 insideConvexHullDir = StereographicProjection::toXyz(insideConvexHullPos);

        Spectrum onConvexHullSp = getSpectrum(onConvexHullDir);
        Spectrum insideConvexHullSp = getSpectrum(insideConvexHullDir);

        if (onConvexHullSp.size() == 0 || insideConvexHullSp.size() == 0) {
            continue;
        }

        Vec2::Scalar t = (v - insideConvexHullPos).norm() / espDist;
        Spectrum     onUnitCircleSp = lerp(insideConvexHullSp, onConvexHullSp, t);
        onUnitCircleSp = onUnitCircleSp.cwiseMax(0.0f);

        sampleMap_[v] = onUnitCircleSp;
        dt_.addVertex(v);
        spectra_.push_back(onUnitCircleSp);
    }

    dt_.computeDelaunayTriangles();
}
