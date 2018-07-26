// =================================================================== //
// Copyright (C) 2014-2018 Kimura Ryo                                  //
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
class SampleSet2D;
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

/*! \brief Fills the samples at the incoming polar angle of 0 degrees using rotational symmetry. */
void fillSpectraAtInThetaOf0(Brdf* brdf);

/*!
 * \brief Rotates a BRDF using an outgoing azimuthal angle.
 * \return A new BRDF with rotated angles.
 */
SphericalCoordinatesBrdf* rotateOutPhi(const SphericalCoordinatesBrdf&  brdf,
                                       float                            rotationAngle);

/*! \brief Fixes the energy conservation of the BRDF with each incoming direction. */
void fixEnergyConservation(SpecularCoordinatesBrdf* brdf);

/*!
 * \brief Fixes the energy conservation of the BRDF with each incoming direction.
 *
 * Reflected energy of a BRDF is reduced if the sum of reflectances of a BRDF and
 * specular reflectance exceed one.
 */
void fixEnergyConservation(SpecularCoordinatesBrdf* brdf,
                           const SampleSet2D&       specularReflectances);

/*!
 * \brief Fixes the energy conservation of the BSDF with each incoming direction.
 *
 * Reflected and transmitted energy of a BRDF and BTDF are reduced if the sum of
 * reflectances and transmittances exceed one.
 */
void fixEnergyConservation(SpecularCoordinatesBrdf* brdf,
                           SpecularCoordinatesBrdf* btdf);

/*! \brief Fills the back side of a BRDF. */
void fillBackSide(SpecularCoordinatesBrdf* brdf);

/*! \brief Removes specular peaks and interpolates using surrounding values. */
void removeSpecularValues(SpecularCoordinatesBrdf* brdf, float maxSpecularTheta);

/*!
 * \brief Insert a BRDF in a base BRDF along incoming azimuthal angle.
 *
 * Both BRDFs must have the same color model, wavelengths, and angles without incoming azimuth.
 *
 * \return A new BRDF with an inserted angle.
 */
Brdf* insertBrdfAlongInPhi(const SphericalCoordinatesBrdf&  baseBrdf,
                           const SphericalCoordinatesBrdf&  insertedBrdf,
                           float                            inPhi);

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
                                         float                          maxSpecularTheta);

/*!
 * \brief Copies spectra from the azimuthal angle of 0 degrees to 360 degrees.
 *
 * Fixes values at the end of azimuthal angles if the start and end points have the same angle.
 */
void copySpectraFromPhiOfZeroTo2PI(SampleSet* samples);

/*! \brief Fills spectra of samples if incoming polar angle is 90 degrees. */
bool fillSpectraAtInThetaOf90(Brdf* brdf, Spectrum::Scalar value = 0.0);

/*! \brief Converts the color model from CIE-XYZ to sRGB. */
void xyzToSrgb(SampleSet* samples);

/*! \brief Fills spectra of samples with a value. */
void fillSpectra(SampleSet* samples, Spectrum::Scalar value);

/*! \brief Fills spectra with a value. */
void fillSpectra(SpectrumList& spectra, Spectrum::Scalar value);

/*! \brief Subtracts a BRDF from another. */
bool subtract(const Brdf& src0, const Brdf& src1, Brdf* dest);

/*! \brief Multiplies spectra of samples by a value. */
void multiplySpectra(SampleSet* samples, Spectrum::Scalar value);

/*! \brief Fixes negative values of spectra to 0. */
void fixNegativeSpectra(SampleSet* samples);

/*! \brief Fixes negative values of spectra to 0. */
void fixNegativeSpectra(SpectrumList& spectra);

} // namespace lb

#endif // LIBBSDF_PROCESSOR_H
