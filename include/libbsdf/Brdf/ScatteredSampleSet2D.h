// =================================================================== //
// Copyright (C) 2022 Kimura Ryo                                       //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#ifndef LIBBSDF_SCATTERED_SAMPLE_SET_2D_H
#define LIBBSDF_SCATTERED_SAMPLE_SET_2D_H

#include <map>

#include <libbsdf/Common/DelaunayTriangulation.h>
#include <libbsdf/Common/Global.h>
#include <libbsdf/Common/StereographicProjection.h>
#include <libbsdf/Common/Utility.h>
#include <libbsdf/Common/Vector.h>

namespace lb {

/*!
 * \class   ScatteredSampleSet2D
 * \brief   The ScatteredSampleSet2D class provides functions for 2D scattered sample points.
 *
 * The data structure consists of pairs of a direction and a spectrum.
 */
class ScatteredSampleSet2D
{
public:
    using SampleMap = std::map<Vec2, Spectrum, CompareVec2>;

    virtual ~ScatteredSampleSet2D();

    /*! Constructs data from added samples. */
    bool constructData(bool extrapolationOfSamplesOnUnitSphere = false);

    /*!
     * Gets the interpolated spectrum at a direction. constructData() must have been called first.
     * An empty spectrum is returned if the direction is outside the convex hull of scattered points.
     */
    Spectrum getSpectrum(const Vec3& dir) const;

    /*! Gets sample points. */
    SampleMap& getSampleMap();

    /*! Gets sample points. */
    const SampleMap& getSampleMap() const;

    /*! Gets Delaunay triangulation of 2D sample points. */
    const DelaunayTriangulation& getDelaunayTriangulation() const;

    /*! Adds a pair of a direction and a spectrum to the list of sample points. */
    void addSample(const Vec3& dir, const Spectrum& sp);

private:
    /*! Computes the intersection of the unit circle with a line starting from inside the circle. */
    Vec2 computeIntersectionOnUnitCircle(const Vec2& pos, const Vec2& dir);

    /*! Adds extrapolated sample points on a unit sphere using preconstructed delaunay triangles. */
    void addExtrapolatedSamplesOnUnitSphere();

    /*!
     * Sample points containing a direction and a spectrum.
     * Directions are transformed by stereographic projection for Delaunay triangulation on a unit sphere.
     */
    SampleMap sampleMap_;

    DelaunayTriangulation dt_; /*!< Delaunay triangulation of 2D sample points. */

    std::vector<Spectrum> spectra_; /*!< The list of spectrum for each direction. */
};

inline ScatteredSampleSet2D::SampleMap& ScatteredSampleSet2D::getSampleMap() { return sampleMap_; }

inline const ScatteredSampleSet2D::SampleMap& ScatteredSampleSet2D::getSampleMap() const
{
    return sampleMap_;
}

inline const DelaunayTriangulation& ScatteredSampleSet2D::getDelaunayTriangulation() const
{
    return dt_;
}

inline void ScatteredSampleSet2D::addSample(const Vec3& dir, const Spectrum& sp)
{
    Vec2 xy = StereographicProjection::fromXyz(dir);
    sampleMap_[xy] = sp;
}

} // namespace lb

#endif // LIBBSDF_SCATTERED_SAMPLE_SET_2D_H
