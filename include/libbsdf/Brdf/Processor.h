// =================================================================== //
// Copyright (C) 2014-2015 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

/*!
 * \file    Processor.h
 * \brief   The Processor.h header file includes the functions to process BRDFs.
 */

#ifndef LIBBSDF_PROCESSOR_H
#define LIBBSDF_PROCESSOR_H

#include <libbsdf/Brdf/Brdf.h>
#include <libbsdf/Brdf/SpecularCoordinatesBrdf.h>
#include <libbsdf/Brdf/SphericalCoordinatesBrdf.h>

namespace lb {

/*!
 * \brief Initializes all spectra of a BRDF using another BRDF.
 *
 * BRDFs must have the same wavelengths.
 */
bool initializeSpectra(const Brdf& baseBrdf, Brdf* brdf);
    
/*!
 * \brief Divides a BRDF by the cosine of the outgoing polar angle.
 *
 * This function converts from a CCBRDF to a BRDF.
 */
void divideByCosineOutTheta(Brdf* brdf);

/*! Fills omitted data using a plane symmetry. */
SphericalCoordinatesBrdf* fillSymmetricBrdf(SphericalCoordinatesBrdf* brdf);

/*! Rotates a BRDF using an outgoing azimuthal angle. */
SphericalCoordinatesBrdf* rotateOutPhi(const SphericalCoordinatesBrdf&  brdf,
                                       float                            rotationAngle);

/*! Fixes the energy conservation of the BRDF with each incoming direction. */
void fixEnergyConservation(SpecularCoordinatesBrdf* brdf);

/*! Copies spectra from the azimuthal angle of 0 to 2PI. */
void copySpectraFromPhiOfZeroTo2PI(Brdf* brdf);

/*! Converts the color model from CIE-XYZ to sRGB. */
void xyzToSrgb(SampleSet* samples);

/*! Fills spectra of samples with a value. */
void fillSpectra(SampleSet* samples, float value);

/*! Multiplies spectra of samples by a value. */
void multiplySpectra(SampleSet* samples, float value);

/*! Clamps negative values of spectra. */
void clampNegativeSpectra(SampleSet* samples);

} // namespace lb

#endif // LIBBSDF_PROCESSOR_H
