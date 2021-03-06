// =================================================================== //
// Copyright (C) 2018-2020 Kimura Ryo                                  //
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

/*! \brief Computes a reflectance of BRDF with a spherical coordinate system. */
Spectrum computeReflectance(const SphericalCoordinatesBrdf& brdf, int inThIndex, int inPhIndex);

/*! \brief Computes a reflectance of BRDF with a specular coordinate system. */
Spectrum computeReflectance(const SpecularCoordinatesBrdf& brdf, int inThIndex, int inPhIndex);

/*!
 * \brief Computes a reflectance of BRDF at an incoming direction.
 *
 *  Partial reflectances are calculated in a specular coordinate system and integrated.
 *
 * \param numThetaDivisions The division number of specular polar angles.
 * \param numPhiDivisions   The division number of specular azimuthal angles.
 */
Spectrum computeReflectance(const Brdf& brdf,
                            const Vec3& inDir,
                            int         numThetaDivisions = 90,
                            int         numPhiDivisions = 72);

/*! \brief Computes reflectances at each incoming direction. */
SampleSet2D* computeReflectances(const SpecularCoordinatesBrdf& brdf);

/*! \brief Computes bihemispherical reflectance (white sky albedo). */
Spectrum computeBihemisphericalReflectance(const Brdf&  brdf,
                                           int          numInThetaDivisions = 9,
                                           int          numInPhiDivisions = 36);

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
 * \brief Computes the bilateral symmetry of BRDF.
 * 
 * The returned spectrum is the bihemispherical reflectance (white sky albedo) of
 * the absolute difference between the original and reversed BRDF.
 *
 * \param numInThetaDivisions   The division number of incoming polar angles to get the bihemispherical reflectance.
 * \param numInPhiDivisions     The division number of incoming azimuthal angles to get the bihemispherical reflectance.
 */
Spectrum computeBilateralSymmetry(const Brdf&   brdf,
                                  int           numInThetaDivisions = 9,
                                  int           numInPhiDivisions = 36);

/*!
 * \brief Computes the reciprocity of BRDF.
 * 
 * The returned spectrum is the bihemispherical reflectance (white sky albedo) of
 * the absolute difference between the original and reversed BRDF.
 *
 * \param numInThetaDivisions   The division number of incoming polar angles to get the bihemispherical reflectance.
 * \param numInPhiDivisions     The division number of incoming azimuthal angles to get the bihemispherical reflectance.
 */
Spectrum computeReciprocity(const Brdf& brdf,
                            int         numInThetaDivisions = 9,
                            int         numInPhiDivisions = 36);

/*!
 * \brief Finds thresholds to separate the diffuse component from a BRDF.
 * \param maxTheta Maxmum incoming and outgoing polar angle to define the range of search.
 */
Spectrum findDiffuseThresholds(const Brdf&      brdf,
                               const double&    maxTheta = PI_2_F);

/*! \brief Returns true if a coordinate system has the angles of an incoming direction. */
bool isInDirDependentCoordinateSystem(const Brdf& brdf);

} // namespace lb

#endif // LIBBSDF_ANALYZER_H
