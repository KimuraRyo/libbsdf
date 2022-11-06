// =================================================================== //
// Copyright (C) 2022 Kimura Ryo                                       //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#include <libbsdf/Common/DelaunayTriangulation.h>

#include <libbsdf/Common/Log.h>

#include <src/ThirdParty/delaunator-cpp/delaunator.hpp>

using namespace lb;

DelaunayTriangulation::DelaunayTriangulation()
{
}

void DelaunayTriangulation::initialize(size_t numVertices)
{
    coords_.clear();
    coords_.reserve(numVertices * 2);

    convexHull_.clear();
}

bool DelaunayTriangulation::computeDelaunayTriangles()
{
    if (coords_.empty()) {
        lbWarn << "[DelaunayTriangulation::computeDelaunayTriangles] Coordinates are empty.";
        return false;
    }

    delaunator::Delaunator d(coords_);

    // Set vertices on convex hull.
    size_t h = d.hull_start;
    while (true) {
        convexHull_.emplace_back(Vec2(d.coords[2 * h], d.coords[2 * h + 1]));
        h = d.hull_next[h];

        if (h == d.hull_start) break;
    }

    triangles_ = std::move(d.triangles);

    return true;
}

std::set<size_t> DelaunayTriangulation::getSurroundingPoints(const double& x, const double& y)
{
    std::set<size_t> indices;

    // Handling of duplicate coordinates
    for (size_t i = 0; i < coords_.size(); i += 2) {
        if (coords_.at(i) == x && coords_.at(i + 1) == y) {
            indices.insert(i);
            indices.insert(i + 1);
            return indices;
        }
    }

    coords_.push_back(x);
    coords_.push_back(y);

    delaunator::Delaunator d(coords_);

    size_t pointIndex = coords_.size() - 2;
    const std::vector<size_t>& trgs = d.triangles;
    for (size_t i = 0; i < trgs.size(); i += 3) {
        size_t index0 = trgs.at(i) * 2;
        size_t index1 = trgs.at(i + 1) * 2;
        size_t index2 = trgs.at(i + 2) * 2;

        if (index0 == pointIndex) {
            indices.insert(index1);
            indices.insert(index2);
        }
        else if (index1 == pointIndex) {
            indices.insert(index2);
            indices.insert(index0);
        }
        else if (index2 == pointIndex) {
            indices.insert(index0);
            indices.insert(index1);
        }
    }

    coords_.pop_back();
    coords_.pop_back();

    return indices;
}
