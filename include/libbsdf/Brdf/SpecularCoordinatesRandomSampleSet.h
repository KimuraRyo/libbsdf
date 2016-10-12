// =================================================================== //
// Copyright (C) 2016 Kimura Ryo                                       //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#ifndef LIBBSDF_SPECULAR_COORDINATES_RANDOM_SAMPLE_SET_H
#define LIBBSDF_SPECULAR_COORDINATES_RANDOM_SAMPLE_SET_H

#include <libbsdf/Brdf/RandomSampleSet.h>
#include <libbsdf/Brdf/SpecularCoordinatesBrdf.h>

namespace lb {

/*!
 * \class   SpecularCoordinatesRandomSampleSet
 * \brief   The SpecularCoordinatesRandomSampleSet class provides the BRDF of a specular coordinate system using random sample points.
 */
class SpecularCoordinatesRandomSampleSet : public RandomSampleSet<SpecularCoordinateSystem>
{
public:
    /*! Sets up a BRDF using random sample points. */
    void setupBrdf(SpecularCoordinatesBrdf* brdf,
                   float                    weight0 = 1.0f,
                   float                    weight1 = 1.0f,
                   float                    weight2 = 1.0f,
                   float                    weight3 = 1.0f);
};

} // namespace lb

#endif // LIBBSDF_SPECULAR_COORDINATES_RANDOM_SAMPLE_SET_H
