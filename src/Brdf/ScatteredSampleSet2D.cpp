// =================================================================== //
// Copyright (C) 2022-2023 Kimura Ryo                                  //
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

    if (extrapolationOfSamplesOnUnitSphere) {
        addExtrapolatedSamplesOnUnitSphere();
    }

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

    return Spectrum();
}

void ScatteredSampleSet2D::addInterpolatedSamples()
{
    const std::vector<size_t>& halfedges = dt_.getHalfEdges();

    std::set<size_t> processedEdges;
    for (size_t edgeIndex = 0; edgeIndex < halfedges.size(); ++edgeIndex) {
        size_t oppositeEdgeIndex = halfedges.at(edgeIndex);

        if (oppositeEdgeIndex == DelaunayTriangulation::INVALID_INDEX){
            size_t edgeVertIndex = dt_.getTriangles().at(edgeIndex);
            size_t oppositeEdgeVertIndex = dt_.getTriangles().at(dt_.getNextHalfEdge(edgeIndex));

            Vec3 vert0(StereographicProjection::toXyz(dt_.getVertex(edgeVertIndex)));
            Vec3 vert1(StereographicProjection::toXyz(dt_.getVertex(oppositeEdgeVertIndex)));

            Vec3 vert = (vert0 + vert1).normalized();
            Vec2 v = StereographicProjection::fromXyz(vert);

            const Spectrum& sp0 = spectra_.at(edgeVertIndex);
            const Spectrum& sp1 = spectra_.at(oppositeEdgeVertIndex);

            Spectrum sp = (sp0 + sp1) * 0.5f;

            sampleMap_[v] = sp;
            dt_.addVertex(v);
            spectra_.push_back(sp);
        }
        else {
            auto ret = processedEdges.insert(edgeIndex);
            if (!ret.second) continue;

            ret = processedEdges.insert(oppositeEdgeIndex);
            if (!ret.second) continue;

            size_t edgeVertIndex = dt_.getTriangles().at(edgeIndex);
            size_t oppositeEdgeVertIndex = dt_.getTriangles().at(oppositeEdgeIndex);

            // The triangle vertex not included in the half-edge
            size_t vert2Index = dt_.getPrevHalfEdge(edgeVertIndex);

            // The opposite triangle vertex not included in the half-edge
            size_t vert3Index = dt_.getPrevHalfEdge(oppositeEdgeVertIndex);

            Vec3 vert0(StereographicProjection::toXyz(dt_.getVertex(edgeVertIndex)));
            Vec3 vert1(StereographicProjection::toXyz(dt_.getVertex(oppositeEdgeVertIndex)));
            Vec3 vert2(StereographicProjection::toXyz(dt_.getVertex(vert2Index)));
            Vec3 vert3(StereographicProjection::toXyz(dt_.getVertex(vert3Index)));

            Vec3 vert = (vert0 + vert1).normalized();
            Vec2 v = StereographicProjection::fromXyz(vert);

            const Spectrum& sp0 = spectra_.at(edgeVertIndex);
            const Spectrum& sp1 = spectra_.at(oppositeEdgeVertIndex);
            const Spectrum& sp2 = spectra_.at(vert2Index);
            const Spectrum& sp3 = spectra_.at(vert3Index);

            double distance0 = std::acos(vert.dot(vert0));
            double distance1 = std::acos(vert.dot(vert1));
            double distance2 = std::acos(vert.dot(vert2));
            double distance3 = std::acos(vert.dot(vert3));

            double weight0 = 1.0 / (distance0 * distance0);
            double weight1 = 1.0 / (distance1 * distance1);
            double weight2 = 1.0 / (distance2 * distance2);
            double weight3 = 1.0 / (distance3 * distance3);

            // Compute a spectrum with inverse distance weighting (IDW).
            Spectrum sumSp = weight0 * sp0 + weight1 * sp1 + weight2 * sp2 + weight3 * sp3;
            double   sumWeight = weight0 + weight1 + weight2 + weight3;
            Spectrum sp = sumSp / sumWeight;

            sampleMap_[v] = sp;
            dt_.addVertex(v);
            spectra_.push_back(sp);
        }
    }

    dt_.computeDelaunayTriangles();
}

Vec2 ScatteredSampleSet2D::computeIntersectionOnCircle(const Vec2&   pos,
                                                       const Vec2&   dir,
                                                       const double& radius)
{
    // https://bluebill.net/circle_ray_intersection/

    Vec2         u = -pos;
    Vec2         u1 = u.dot(dir) * dir;
    Vec2         u2 = u - u1;
    Vec2::Scalar d = u2.norm();
    Vec2::Scalar m = std::sqrt(radius * radius - d * d);

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
    const std::vector<size_t>& convexHull = dt_.getConvexHull();

    for (size_t i = 0; i < convexHull.size(); ++i) {
        size_t backI = convexHull.size() - 1;
        size_t prevI = (i != 0) ? i - 1 : backI;
        size_t nextI = (i != backI) ? i + 1 : 0;

        size_t h = convexHull.at(i);
        size_t prevH = convexHull.at(prevI);
        size_t nextH = convexHull.at(nextI);

        Vec2 pos = dt_.getVertex(h);

        Vec2 prevEdgeDir = (pos - dt_.getVertex(prevH)).normalized();
        Vec2 nextEdgeDir = (pos - dt_.getVertex(nextH)).normalized();

        if (prevEdgeDir.dot(nextEdgeDir) < -0.99999) continue;

        Vec2 dir = (prevEdgeDir + nextEdgeDir).normalized();

        // Slightly increase the radius of the stereographic projected circle
        // to get a spectrum at a polar angle of 90 degrees.
        constexpr double radius = 1.1;
        Vec2 posOnCircle = computeIntersectionOnCircle(pos, dir, radius);

        constexpr Vec2::Scalar espDist = EPSILON_F;

        Vec2 posInsideConvexHull = pos - espDist * dir;

        Vec3 dirOnConvexHull = StereographicProjection::toXyz(pos);
        Vec3 dirInsideConvexHull = StereographicProjection::toXyz(posInsideConvexHull);

        Spectrum spOnConvexHull = spectra_.at(h);
        Spectrum spInsideConvexHull = getSpectrum(dirInsideConvexHull);

        if (spOnConvexHull.size() == 0 || spInsideConvexHull.size() == 0) {
            lbWarn << "[ScatteredSampleSet2D::addExtrapolatedSamplesOnUnitSphere] Failed to extrapolate.";
            continue;
        }

        Vec2::Scalar t = (posOnCircle - posInsideConvexHull).norm() / espDist;
        Spectrum     spOnCircle = lerp(spInsideConvexHull, spOnConvexHull, t);

        sampleMap_[posOnCircle] = spOnCircle;
        dt_.addVertex(posOnCircle);
        spectra_.push_back(spOnCircle);
    }

    dt_.computeDelaunayTriangles();
}
