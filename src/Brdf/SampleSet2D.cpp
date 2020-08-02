// =================================================================== //
// Copyright (C) 2014-2020 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#include <libbsdf/Brdf/SampleSet2D.h>

#include <libbsdf/Brdf/Sampler.h>

using namespace lb;

SampleSet2D::SampleSet2D(int        numTheta,
                         int        numPhi,
                         ColorModel colorModel,
                         int        numWavelengths,
                         bool       equalIntervalAngles)
                         : equalIntervalTheta_(false),
                           equalIntervalPhi_(false),
                           sourceType_(UNKNOWN_SOURCE)
{
    lbTrace << "[SampleSet2D::SampleSet2D]";

    assert(numTheta > 0 && numPhi > 0);

    colorModel_ = colorModel;

    size_t numSamples = numTheta * numPhi;
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

    for (size_t i = 0; i < numSamples; ++i) {
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

SampleSet2D::~SampleSet2D()
{
    lbTrace << "[SampleSet2D::~SampleSet2D]";
}

bool SampleSet2D::validate(bool verbose) const
{
    bool valid = true;
    bool spectraValid = true;

    // Spectra
    for (int i0 = 0; i0 < thetaAngles_.size(); ++i0) {
        if (!spectraValid && !verbose) break;
    for (int i1 = 0; i1 < phiAngles_.size(); ++i1) {
        const Spectrum& sp = getSpectrum(i0, i1);

        if (!sp.allFinite()) {
            spectraValid = false;

            if (sp.hasNaN()) {
                lbWarn
                    << "[SampleSet2D::validate] The spectrum contains NaN value(s) at ("
                    << i0 << ", " << i1 << "):\n\t"
                    << sp.format(LB_EIGEN_IO_FMT);
            }
            else {
                lbWarn
                    << "[SampleSet2D::validate] The spectrum contains +/-INF value(s) at ("
                    << i0 << ", " << i1 << "):\n\t"
                    << sp.format(LB_EIGEN_IO_FMT);
            }

            if (!verbose) break;
        }
    }}

    if (spectraValid) {
        lbInfo << "[SampleSet2D::validate] Spectra are valid.";
    }
    else {
        valid = false;
        lbWarn << "[SampleSet2D::validate] Invalid spectra are found.";
    }

    // Angle arrays
    if (thetaAngles_.allFinite()) {
        lbInfo << "[SampleSet2D::validate] The array of angle0 is valid.";
    }
    else {
        valid = false;
        lbWarn << "[SampleSet2D::validate] The invalid angle(s) in angles0 is found.";
    }

    if (phiAngles_.allFinite()) {
        lbInfo << "[SampleSet2D::validate] The array of angle1 is valid.";
    }
    else {
        valid = false;
        lbWarn << "[SampleSet2D::validate] The invalid angle(s) in angles1 is found.";
    }

    // Angle attributes
    if (equalIntervalTheta_ && !array_util::isEqualInterval(thetaAngles_)) {
        valid = false;
        lbWarn << "[SampleSet2D::validate] equalIntervalTheta_ attribute is not updated.";
    }

    if (equalIntervalPhi_ && !array_util::isEqualInterval(phiAngles_)) {
        valid = false;
        lbWarn << "[SampleSet2D::validate] equalIntervalPhi_ attribute is not updated.";
    }

    // Wavelengths
    if (wavelengths_.allFinite()) {
        if (wavelengths_.minCoeff() < 0.0f) {
            lbWarn
                << "[SampleSet2D::validate] The negative wavelength(s) is found:\n\t"
                << wavelengths_.format(LB_EIGEN_IO_FMT);
        }
        else {
            lbInfo << "[SampleSet2D::validate] Wavelengths are valid.";
        }
    }
    else {
        valid = false;
        lbWarn
            << "[SampleSet2D::validate] The invalid wavelength(s) is found:\n\t"
            << wavelengths_.format(LB_EIGEN_IO_FMT);
    }

    return valid;
}

Spectrum SampleSet2D::getSpectrum(const Vec3& inDir) const
{
    return Sampler::getSpectrum<LinearInterpolator>(*this, inDir);
}

void SampleSet2D::updateAngleAttributes()
{
    equalIntervalTheta_ = array_util::isEqualInterval(thetaAngles_);
    equalIntervalPhi_   = array_util::isEqualInterval(phiAngles_);

    lbInfo
        << "[SampleSet2D::updateAngleAttributes] equalIntervalTheta_: "
        << equalIntervalTheta_;
    lbInfo
        << "[SampleSet2D::updateAngleAttributes] equalIntervalPhi_: "
        << equalIntervalPhi_;
}

void SampleSet2D::resizeAngles(int numTheta, int numPhi)
{
    assert(numTheta > 0 && numPhi > 0);

    size_t numSamples = numTheta * numPhi;
    spectra_.resize(numSamples);

    thetaAngles_.resize(numTheta);
    phiAngles_.resize(numPhi);
}

void SampleSet2D::resizeWavelengths(int numWavelengths)
{
    assert(numWavelengths > 0);

    size_t numSamples = thetaAngles_.size() * phiAngles_.size();

    for (size_t i = 0; i < numSamples; ++i) {
        Spectrum sp;
        sp.resize(numWavelengths);
        spectra_.at(i) = sp;
    }

    wavelengths_.resize(numWavelengths);
}

void SampleSet2D::clampAngles()
{
    thetaAngles_ = thetaAngles_.cwiseMax(0.0);
    phiAngles_   = phiAngles_.cwiseMax(0.0);

    thetaAngles_ = thetaAngles_.cwiseMin(SphericalCoordinateSystem::MAX_ANGLE0);
    phiAngles_   = phiAngles_.cwiseMin(SphericalCoordinateSystem::MAX_ANGLE1);
}
