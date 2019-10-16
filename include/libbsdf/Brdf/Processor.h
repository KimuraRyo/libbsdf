// =================================================================== //
// Copyright (C) 2014-2019 Kimura Ryo                                  //
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

#include <functional>

#include <libbsdf/Brdf/SampleSet2D.h>
#include <libbsdf/Brdf/SpecularCoordinatesBrdf.h>
#include <libbsdf/Brdf/SphericalCoordinatesBrdf.h>

namespace lb {

/*!
 * \brief Edits diffuse and glossy components of a BRDF.
 * \param diffuseThresholds Thresholds to separate the diffuse component from a BRDF.
 * \param glossyIntensity   Rate of change of a glossy intensity. If 1.0 is specified, it is not changed.
 * \param glossyShininess   Ratio of the thinness of a glossy lobe. If 1.0 is specified, shininess is not changed.
 * \param diffuseIntensity  Rate of change of a diffuse intensity. If 1.0 is specified, it is not changed.
 */
void editComponents(const Brdf&         origBrdf,
                    Brdf*               brdf,
                    const Spectrum&     diffuseThresholds,
                    Spectrum::Scalar    glossyIntensity,
                    Spectrum::Scalar    glossyShininess,
                    Spectrum::Scalar    diffuseIntensity);

/*!
 * \brief Edits diffuse and glossy components of a BRDF value at a sample point.
 * \param diffuseThresholds Thresholds to separate the diffuse component from a BRDF.
 * \param glossyIntensity   Rate of change of a glossy intensity. If 1.0 is specified, it is not changed.
 * \param glossyShininess   Ratio of the thinness of a glossy lobe. If 1.0 is specified, shininess is not changed.
 * \param diffuseIntensity  Rate of change of a diffuse intensity. If 1.0 is specified, it is not changed.
 */
void editComponents(int                 i0,
                    int                 i1,
                    int                 i2,
                    int                 i3,
                    const Brdf&         origBrdf,
                    Brdf*               brdf,
                    const Spectrum&     diffuseThresholds,
                    Spectrum::Scalar    glossyIntensity,
                    Spectrum::Scalar    glossyShininess,
                    Spectrum::Scalar    diffuseIntensity);

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

/*!
 * \brief Fills the back side of a BRDF.
 *
 * Samples in the back side of a surface are filled with copies of boundary samples in the front side.
 */
void fillBackSide(SpecularCoordinatesBrdf* brdf);

/*! \brief Averages samples for overlapping angles. */
void equalizeOverlappingSamples(SpecularCoordinatesBrdf* brdf);

/*!
 * \brief Removes specular peaks and interpolates removed areas using surrounding values.
 *
 * \param maxSpecularTheta  The maximum angle of removed areas around specular directions.
 */
void removeSpecularValues(SpecularCoordinatesBrdf* brdf, float maxSpecularTheta);

/*!
 * \brief Inserts a BRDF in a base BRDF along incoming azimuthal angle.
 *
 * Both BRDFs must have the same color model, wavelengths, and angles without incoming azimuth.
 *
 * \return A new BRDF with a specular coordinate system.
 */
SpecularCoordinatesBrdf* insertBrdfAlongInPhi(const SpecularCoordinatesBrdf&   baseBrdf,
                                              const SpecularCoordinatesBrdf&   insertedBrdf,
                                              float                            inPhi);

/*!
 * \brief Inserts a BRDF in a base BRDF along incoming azimuthal angle.
 *
 * Both BRDFs must have the same color model, wavelengths, and angles without incoming azimuth.
 *
 * \return A new BRDF with a spherical coordinate system.
 */
SphericalCoordinatesBrdf* insertBrdfAlongInPhi(const SphericalCoordinatesBrdf&  baseBrdf,
                                               const SphericalCoordinatesBrdf&  insertedBrdf,
                                               float                            inPhi);

/*!
 * \brief Recalculates a BRDF with linearly extrapolated reflectances.
 * \param incomingTheta Minimum incoming polar angle of extrapolated samples.
 *
 * When a BRDF is extrapolated, it is separated into glossy and diffuse components to improve calculation results.
 */
void extrapolateSamplesWithReflectances(SpecularCoordinatesBrdf* brdf, float incomingTheta);

/*!
 * \brief Copies spectra from the azimuthal angle of 0 degrees to 360 degrees.
 *
 * Fixes values at the end of azimuthal angles if the start and end points have the same angle.
 */
void copySpectraFromPhiOf0To360(SampleSet* samples);

/*! \brief Fills spectra of samples if incoming polar angle is 90 degrees. */
bool fillSpectraAtInThetaOf90(Brdf* brdf, Spectrum::Scalar value = 0.0);

/*! \brief Converts the color model from CIE-XYZ to sRGB. */
void xyzToSrgb(SampleSet* samples);

/*! \brief Fills spectra of samples with a value. */
void fillSpectra(SampleSet* samples, Spectrum::Scalar value);

/*! \brief Fills spectra with a value. */
void fillSpectra(SpectrumList& spectra, Spectrum::Scalar value);

/*! \brief Computes spectra of a BRDF with two BRDFs. */
bool compute(const Brdf& src0, const Brdf& src1, Brdf* dest,
             std::function<Spectrum(const Spectrum&, const Spectrum&)> manipulator);

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
