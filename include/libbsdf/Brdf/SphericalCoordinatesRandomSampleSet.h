// =================================================================== //
// Copyright (C) 2016 Kimura Ryo                                       //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#ifndef LIBBSDF_SPHERICAL_COORDINATES_RANDOM_SAMPLE_SET_H
#define LIBBSDF_SPHERICAL_COORDINATES_RANDOM_SAMPLE_SET_H

#include <libbsdf/Brdf/RandomSampleSet.h>
#include <libbsdf/Brdf/SphericalCoordinatesBrdf.h>

namespace lb {

/*!
 * \class   SphericalCoordinatesRandomSampleSet
 * \brief   The SphericalCoordinatesRandomSampleSet class provides the BRDF of a spherical coordinate system using random sample points.
 */
class SphericalCoordinatesRandomSampleSet : public RandomSampleSet<SphericalCoordinateSystem>
{
public:
    /*! Sets up a BRDF using random sample points. */
    void setupBrdf(SphericalCoordinatesBrdf* brdf);
};

} // namespace lb

#endif // LIBBSDF_SPHERICAL_COORDINATES_RANDOM_SAMPLE_SET_H
