// =================================================================== //
// Copyright (C) 2014-2015 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#include <libbsdf/Brdf/SampleSet2D.h>

#include <libbsdf/Brdf/Sampler.h>

#include <iostream>

using namespace lb;

SampleSet2D::SampleSet2D(int        numTheta,
                         int        numPhi,
                         ColorModel colorModel,
                         int        numWavelengths,
                         bool       equalIntervalAngles)
                         : equalIntervalTheta_(false),
                           equalIntervalPhi_(false)
{
    assert(numTheta > 0 && numPhi > 0);

    numTheta_ = numTheta;
    numPhi_ = numPhi;

    colorModel_ = colorModel;

    int numSamples = numTheta * numPhi;
    spectra_.resize(numSamples);

    if (equalIntervalAngles) {
        thetaAngles_ = Arrayf::LinSpaced(numTheta, 0.0, SphericalCoordinateSystem::MAX_ANGLE0);
        phiAngles_   = Arrayf::LinSpaced(numPhi,   0.0, SphericalCoordinateSystem::MAX_ANGLE1);
        updateAngleAttributes();
    }
    else {
        thetaAngles_.resize(numTheta);
        phiAngles_.resize(numPhi);
    }

    if (colorModel == MONOCHROMATIC_MODEL) {
        numWavelengths = 1;
    }
    else if (colorModel != SPECTRAL_MODEL) {
        numWavelengths = 3;
    }

    assert(numWavelengths > 0);

    for (int i = 0; i < numSamples; ++i) {
        Spectrum spectrum;
        spectrum.resize(numWavelengths);
        spectra_.at(i) = spectrum;
    }

    wavelengths_.resize(numWavelengths);

    if (colorModel == MONOCHROMATIC_MODEL ||
        colorModel != SPECTRAL_MODEL) {
        wavelengths_ = Arrayf::Zero(numWavelengths);
    }
}

Spectrum SampleSet2D::getSpectrum(const Vec3& inDir) const
{
    Spectrum sp;
    Sampler::getSpectrum<LinearInterpolator>(*this, inDir, &sp);
    return sp;
}

void SampleSet2D::updateAngleAttributes()
{
    equalIntervalTheta_ = isEqualInterval(thetaAngles_);
    equalIntervalPhi_   = isEqualInterval(phiAngles_);

    std::cout << "[SampleSet2D::updateAngleAttributes] Theta: " << equalIntervalTheta_ << std::endl;
    std::cout << "[SampleSet2D::updateAngleAttributes] Phi: "   << equalIntervalPhi_   << std::endl;
}

void SampleSet2D::clampAngles()
{
    thetaAngles_ = thetaAngles_.cwiseMax(0.0);
    phiAngles_   = phiAngles_.cwiseMax(0.0);

    thetaAngles_ = thetaAngles_.cwiseMin(SphericalCoordinateSystem::MAX_ANGLE0);
    phiAngles_   = phiAngles_.cwiseMin(SphericalCoordinateSystem::MAX_ANGLE1);
}
