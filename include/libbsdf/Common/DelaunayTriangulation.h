// =================================================================== //
// Copyright (C) 2022 Kimura Ryo                                       //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#ifndef LIBBSDF_DELAUNAY_TRIANGULATION_H
#define LIBBSDF_DELAUNAY_TRIANGULATION_H

#include <cmath>
#include <set>

#include <libbsdf/Common/Global.h>
#include <libbsdf/Common/Vector.h>

namespace lb {

/*!
 * \class   DelaunayTriangulation
 * \brief   The DelaunayTriangulation class provides functions for Delaunay triangulation based on Delaunator.
 *
 * The data structure of Delaunay triangulation is the same as Delaunator.
 * https://mapbox.github.io/delaunator/
 * https://github.com/abellgithub/delaunator-cpp
 */
class DelaunayTriangulation
{
public:
    explicit DelaunayTriangulation();

    /*! Clears and reserves the buffer of vertices. */
    void initialize(size_t numVertices);

    /*! Adds a vertex. */
    void addVertex(const Vec2& position);

    /*! Computes Delaunay triangles and the convex hull. */
    bool computeDelaunayTriangles();

    /*! Gets surrounding points if a point is added. */
    std::set<size_t> getSurroundingPoints(const double& x, const double& y);

    /*! Gets the number of Delaunay triangles. */
    size_t getNumTriangles() const;

    /*! Gets the vertex indices of a triangle. */
    Vec3i getTriangle(const size_t& index) const;

    /*! Gets a vertex. */
    Vec2 getVertex(const size_t& index) const;

    std::vector<double>&       getCoords();
    const std::vector<double>& getCoords() const;
    const std::vector<size_t>& getTriangles() const;

    /*! Gets the vertex positions of the convex hull. */
    const std::vector<Vec2>& getConvexHull() const;

private:
    std::vector<double> coords_;     /*!< 2D coordinates (x0, y0, x1, y1, ...) */
    std::vector<size_t> triangles_;  /*!< Indices to the 'X's of \a coords_ */
    std::vector<Vec2>   convexHull_; /*!< Vertex positions of the convex hull. */
};

inline void DelaunayTriangulation::addVertex(const Vec2& position)
{
    coords_.push_back(position[0]);
    coords_.push_back(position[1]);
}

inline size_t DelaunayTriangulation::getNumTriangles() const { return triangles_.size() / 3; }

inline Vec3i DelaunayTriangulation::getTriangle(const size_t& index) const
{
    using Scalar = Vec3i::Scalar; 
    size_t i = index * 3;
    return Vec3i(static_cast<Scalar>(triangles_.at(i)),
                 static_cast<Scalar>(triangles_.at(i + 1)),
                 static_cast<Scalar>(triangles_.at(i + 2)));
}

inline Vec2 DelaunayTriangulation::getVertex(const size_t& index) const
{
    size_t i = index * 2;
    return Vec2(coords_.at(i), coords_.at(i + 1));
}

inline std::vector<double>&       DelaunayTriangulation::getCoords() { return coords_; }
inline const std::vector<double>& DelaunayTriangulation::getCoords() const { return coords_; }
inline const std::vector<size_t>& DelaunayTriangulation::getTriangles() const { return triangles_; }

inline const std::vector<Vec2>& DelaunayTriangulation::getConvexHull() const { return convexHull_; }

} // namespace lb

#endif // LIBBSDF_DELAUNAY_TRIANGULATION_H
