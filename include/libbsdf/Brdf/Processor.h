// =================================================================== //
// Copyright (C) 2014-2016 Kimura Ryo                                  //
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

#include <libbsdf/Common/Global.h>

namespace lb {

class Brdf;
class SampleSet;
class SpecularCoordinatesBrdf;
class SphericalCoordinatesBrdf;

/*!
 * \brief Divides a BRDF by the cosine of the outgoing polar angle.
 *
 * This function converts from a CCBRDF to a BRDF.
 */
void divideByCosineOutTheta(Brdf* brdf);

/*!
 * \brief Fills omitted data using plane symmetry.
 * \return A new BRDF with appended angles.
 */
SphericalCoordinatesBrdf* fillSymmetricBrdf(SphericalCoordinatesBrdf* brdf);

/*! \brief Fills the samples at the incoming polar angle of 0 using rotational symmetry. */
void fillIncomingPolar0Data(Brdf* brdf);

/*!
 * \brief Rotates a BRDF using an outgoing azimuthal angle.
 * \return A new BRDF with rotated angles.
 */
SphericalCoordinatesBrdf* rotateOutPhi(const SphericalCoordinatesBrdf&  brdf,
                                       float                            rotationAngle);

/*! \brief Fixes the energy conservation of the BRDF with each incoming direction. */
void fixEnergyConservation(SpecularCoordinatesBrdf* brdf);

/*! \brief Copies spectra from the azimuthal angle of 0 to 2PI. */
void copySpectraFromPhiOfZeroTo2PI(Brdf* brdf);

/*! \brief Converts the color model from CIE-XYZ to sRGB. */
void xyzToSrgb(SampleSet* samples);

/*! \brief Fills spectra of samples with a value. */
void fillSpectra(SampleSet* samples, Spectrum::Scalar value);

/*! \brief Fills spectra with a value. */
void fillSpectra(SpectrumList& spectra, Spectrum::Scalar value);

/*! \brief Multiplies spectra of samples by a value. */
void multiplySpectra(SampleSet* samples, Spectrum::Scalar value);

/*! \brief Fixes negative values of spectra to 0. */
void fixNegativeSpectra(SampleSet* samples);

} // namespace lb

#endif // LIBBSDF_PROCESSOR_H
