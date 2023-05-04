// =================================================================== //
// Copyright (C) 2022-2023 Kimura Ryo                                  //
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
    Vec3i getTriangle(size_t index) const;

    /*! Gets the number of vertices. */
    size_t getNumVertices() const;

    /*! Gets a vertex. */
    Vec2 getVertex(size_t index) const;

    std::vector<double>&       getCoords();
    const std::vector<double>& getCoords() const;
    const std::vector<size_t>& getTriangles() const;
    const std::vector<size_t>& getHalfEdges() const;

    /*! Gets the vertex positions of the convex hull. */
    const std::vector<size_t>& getConvexHull() const;

    /*! Gets the triangle index of an edge. */
    size_t getTriangleOfEdge(size_t edgeIndex) const;

    /*! Gets the next half-edge of a half-edge in the same triangle. */
    size_t getNextHalfEdge(size_t edgeIndex) const;

    /*! Gets the previous half-edge of a half-edge in the same triangle. */
    size_t getPrevHalfEdge(size_t edgeIndex) const;

    /*! delaunator::INVALID_INDEX */
    static const size_t INVALID_INDEX;

private:
    std::vector<double> coords_;     /*!< 2D coordinates (x0, y0, x1, y1, ...) */
    std::vector<size_t> triangles_;  /*!< Indices to the 'X's of \a coords_ */
    std::vector<size_t> halfedges_;  /*!< Half-edges. See delaunator.hpp. */
    std::vector<size_t> convexHull_; /*!< Vertex indices of the convex hull. */
};

inline void DelaunayTriangulation::addVertex(const Vec2& position)
{
    coords_.push_back(position[0]);
    coords_.push_back(position[1]);
}

inline size_t DelaunayTriangulation::getNumTriangles() const { return triangles_.size() / 3; }

inline Vec3i DelaunayTriangulation::getTriangle(size_t index) const
{
    using Scalar = Vec3i::Scalar; 
    size_t i = index * 3;
    return Vec3i(static_cast<Scalar>(triangles_.at(i)),
                 static_cast<Scalar>(triangles_.at(i + 1)),
                 static_cast<Scalar>(triangles_.at(i + 2)));
}

inline size_t DelaunayTriangulation::getNumVertices() const
{
    return coords_.size() / 2;
}

inline Vec2 DelaunayTriangulation::getVertex(size_t index) const
{
    size_t i = index * 2;
    return Vec2(coords_.at(i), coords_.at(i + 1));
}

inline std::vector<double>&       DelaunayTriangulation::getCoords() { return coords_; }
inline const std::vector<double>& DelaunayTriangulation::getCoords() const { return coords_; }

inline const std::vector<size_t>& DelaunayTriangulation::getTriangles() const { return triangles_; }
inline const std::vector<size_t>& DelaunayTriangulation::getHalfEdges() const { return halfedges_; }
inline const std::vector<size_t>& DelaunayTriangulation::getConvexHull() const { return convexHull_; }


inline size_t DelaunayTriangulation::getTriangleOfEdge(size_t edgeIndex) const
{
    return static_cast<size_t>(edgeIndex / 3);
}

inline size_t DelaunayTriangulation::getNextHalfEdge(size_t edgeIndex) const
{
    return (edgeIndex % 3 == 2) ? edgeIndex - 2 : edgeIndex + 1;
}

inline size_t DelaunayTriangulation::getPrevHalfEdge(size_t edgeIndex) const
{
    return (edgeIndex % 3 == 0) ? edgeIndex + 2 : edgeIndex - 1;
}

} // namespace lb

#endif // LIBBSDF_DELAUNAY_TRIANGULATION_H
