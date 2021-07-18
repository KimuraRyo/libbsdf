// =================================================================== //
// Copyright (C) 2019-2020 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

/*!
 * \file    Initializer.h
 * \brief   The Initializer.h header file includes the functions to initialize data.
 */

#ifndef LIBBSDF_INITIALIZER_H
#define LIBBSDF_INITIALIZER_H

#include <libbsdf/Brdf/Sampler.h>

namespace lb {

/*!
 * \brief Initializes all spectra of a BRDF using another BRDF.
 *
 * Both BRDFs must have the same color model and wavelengths.
 */
template <typename InterpolatorT>
bool initializeSpectra(const Brdf& baseBrdf, Brdf* brdf);

/*!
 * \brief Initializes all spectra of samples using another samples.
 *
 * Both samples must have the same color model and wavelengths.
 */
template <typename InterpolatorT>
bool initializeSpectra(const SampleSet2D& baseSamples, SampleSet2D* samples);

} // namespace lb

template <typename InterpolatorT>
bool lb::initializeSpectra(const Brdf& baseBrdf, Brdf* brdf)
{
    lbTrace << "[initializeSpectra]";

    const SampleSet* baseSs = baseBrdf.getSampleSet();
    SampleSet* ss = brdf->getSampleSet();

    if (!hasSameColor(*baseSs, *ss)) {
        lbError << "[initializeSpectra] Color models or wavelengths do not match.";
        return false;
    }

    for (int i0 = 0; i0 < ss->getNumAngles0(); ++i0) {
    for (int i1 = 0; i1 < ss->getNumAngles1(); ++i1) {
    for (int i2 = 0; i2 < ss->getNumAngles2(); ++i2) {
        Vec3 inDir, outDir;
        Spectrum sp;
        #pragma omp parallel for private(inDir, outDir, sp)
        for (int i3 = 0; i3 < ss->getNumAngles3(); ++i3) {
            brdf->getInOutDirection(i0, i1, i2, i3, &inDir, &outDir);
            sp = Sampler::getSpectrum<InterpolatorT>(baseBrdf, inDir, outDir);
            ss->setSpectrum(i0, i1, i2, i3, sp);
        }
    }}}

    return true;
}

template <typename InterpolatorT>
bool lb::initializeSpectra(const SampleSet2D& baseSamples, SampleSet2D* samples)
{
    lbTrace << "[initializeSpectra]";

    if (!hasSameColor(baseSamples, *samples)) {
        lbError << "[initializeSpectra] Color models or wavelengths do not match.";
        return false;
    }

    for (int i0 = 0; i0 < samples->getNumTheta(); ++i0) {
    for (int i1 = 0; i1 < samples->getNumPhi();   ++i1) {
        float theta = samples->getTheta(i0);
        float phi   = samples->getPhi(i1);
        Spectrum sp = InterpolatorT::getSpectrum(baseSamples, theta, phi);
        samples->setSpectrum(i0, i1, sp);
    }}

    return true;
}

#endif // LIBBSDF_INITIALIZER_H
