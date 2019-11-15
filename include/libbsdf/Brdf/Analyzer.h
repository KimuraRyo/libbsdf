// =================================================================== //
// Copyright (C) 2018-2019 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

/*!
 * \file    Analyzer.h
 * \brief   The Analyzer.h header file includes the functions to analyze BRDFs.
 */

#ifndef LIBBSDF_ANALYZER_H
#define LIBBSDF_ANALYZER_H

#include <libbsdf/Brdf/SampleSet2D.h>
#include <libbsdf/Brdf/SpecularCoordinatesBrdf.h>
#include <libbsdf/Brdf/SphericalCoordinatesBrdf.h>

namespace lb {

/*!
 * \brief Computes a reflectance of BRDF with a spherical coordinate system.
 */
Spectrum computeReflectance(const SphericalCoordinatesBrdf& brdf, int inThIndex, int inPhIndex);

/*!
 * \brief Computes a reflectance of BRDF with a specular coordinate system.
 */
Spectrum computeReflectance(const SpecularCoordinatesBrdf& brdf, int inThIndex, int inPhIndex);

/*!
 * \brief Computes a reflectance of BRDF at an incoming direction.
 */
Spectrum computeReflectance(const Brdf& brdf, const Vec3& inDir);

/*!
 * \brief Computes reflectances at each incoming direction.
 */
SampleSet2D* computeReflectances(const SpecularCoordinatesBrdf& brdf);

/*!
 * \brief Computes specular reflectances using a standard sample.
 * \param ior   Index of refraction of the standard material. 1.0 is used for transmittance.
 */
SampleSet2D* computeSpecularReflectances(const Brdf&    brdf,
                                         const Brdf&    standardBrdf,
                                         float          ior);

/*!
 * \brief Computes specular reflectances using a standard sample.
 *
 * The maximum spectrum found around a specular direction is used as the specular
 * reflected or transmitted spectrum for each incoming direction. A sample may lean
 * to one side, when it is measured.
 *
 * \param ior               Index of refraction of the standard material. 1.0 is used for transmittance.
 * \param maxSpecularTheta  The maximum angle in radians to find the specular reflected/transmitted spectrum.
 */
SampleSet2D* computeSpecularReflectances(const SpecularCoordinatesBrdf& brdf,
                                         const Brdf&                    standardBrdf,
                                         float                          ior,
                                         float                          maxSpecularTheta = PI_2_F);

/*!
 * \brief Finds thresholds to separate the diffuse component from a BRDF.
 * \param maxTheta Maxmum incoming and outgoing polar angle to define the range of search.
 */
Spectrum findDiffuseThresholds(const Brdf&      brdf,
                               const double&    maxTheta = PI_2_F);

/*! Returns true if a coordinate system has the angles of an incoming direction. */
bool isInDirDependentCoordinateSystem(const Brdf& brdf);

} // namespace lb

#endif // LIBBSDF_ANALYZER_H
