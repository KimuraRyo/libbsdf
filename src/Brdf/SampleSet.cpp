// =================================================================== //
// Copyright (C) 2014-2023 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#include <libbsdf/Brdf/SampleSet.h>

using namespace lb;

SampleSet::SampleSet(int        numAngles0,
                     int        numAngles1,
                     int        numAngles2,
                     int        numAngles3,
                     ColorModel colorModel,
                     int        numWavelengths)
    : equalIntervalAngles0_(false),
      equalIntervalAngles1_(false),
      equalIntervalAngles2_(false),
      equalIntervalAngles3_(false),
      oneSide_(false)
{
    assert(numAngles0 > 0 && numAngles1 > 0 && numAngles2 > 0 && numAngles3 > 0);

    resizeAngles(numAngles0, numAngles1, numAngles2, numAngles3);

    colorModel_ = colorModel;

    if (colorModel == SPECTRAL_MODEL) {
        resizeWavelengths(numWavelengths);
    }
    else if (colorModel == MONOCHROMATIC_MODEL) {
        resizeWavelengths(1);
        wavelengths_ = Arrayf::Zero(1);
    }
    else {
        resizeWavelengths(3);
        wavelengths_ = Arrayf::Zero(3);
    }
}

bool SampleSet::validate(bool verbose) const
{
    bool valid = true;
    bool spectraValid = true;

    // Spectra
    for (int i0 = 0; i0 < angles0_.size(); ++i0) {
        if (!spectraValid && !verbose) break;
    for (int i1 = 0; i1 < angles1_.size(); ++i1) {
        if (!spectraValid && !verbose) break;
    for (int i2 = 0; i2 < angles2_.size(); ++i2) {
        if (!spectraValid && !verbose) break;
    for (int i3 = 0; i3 < angles3_.size(); ++i3) {
        const Spectrum& sp = getSpectrum(i0, i1, i2, i3);

        if (!sp.allFinite()) {
            spectraValid = false;

            if (sp.hasNaN()) {
                lbWarn
                    << "[SampleSet::validate] The spectrum contains NaN value(s) at ("
                    << i0 << ", " << i1 << ", " << i2 << ", " << i3 << "):\n\t"
                    << sp.format(LB_EIGEN_IO_FMT);
            }
            else {
                lbWarn
                    << "[SampleSet::validate] The spectrum contains +/-INF value(s) at ("
                    << i0 << ", " << i1 << ", " << i2 << ", " << i3 << "):\n\t"
                    << sp.format(LB_EIGEN_IO_FMT);
            }

            if (!verbose) break;
        }
    }}}}

    if (spectraValid) {
        lbInfo << "[SampleSet::validate] Spectra are valid.";
    }
    else {
        valid = false;
        lbWarn << "[SampleSet::validate] Invalid spectra are found.";
    }

    // Angle arrays
    if (angles0_.allFinite()) {
        lbInfo << "[SampleSet::validate] The array of angle0 is valid.";
    }
    else {
        valid = false;
        lbWarn << "[SampleSet::validate] The invalid angle(s) in angles0 is found.";
    }

    if (angles1_.allFinite()) {
        lbInfo << "[SampleSet::validate] The array of angle1 is valid.";
    }
    else {
        valid = false;
        lbWarn << "[SampleSet::validate] The invalid angle(s) in angles1 is found.";
    }

    if (angles2_.allFinite()) {
        lbInfo << "[SampleSet::validate] The array of angle2 is valid.";
    }
    else {
        valid = false;
        lbWarn << "[SampleSet::validate] The invalid angle(s) in angles2 is found.";
    }

    if (angles3_.allFinite()) {
        lbInfo << "[SampleSet::validate] The array of angle3 is valid.";
    }
    else {
        valid = false;
        lbWarn << "[SampleSet::validate] The invalid angle(s) in angles3 is found.";
    }

    // Angle attributes
    if (equalIntervalAngles0_ && !array_util::isEqualInterval(angles0_)) {
        valid = false;
        lbWarn << "[SampleSet::validate] equalIntervalAngles0 attribute is not updated.";
    }

    if (equalIntervalAngles1_ && !array_util::isEqualInterval(angles1_)) {
        valid = false;
        lbWarn << "[SampleSet::validate] equalIntervalAngles1 attribute is not updated.";
    }

    if (equalIntervalAngles2_ && !array_util::isEqualInterval(angles2_)) {
        valid = false;
        lbWarn << "[SampleSet::validate] equalIntervalAngles2 attribute is not updated.";
    }

    if (equalIntervalAngles3_ && !array_util::isEqualInterval(angles3_)) {
        valid = false;
        lbWarn << "[SampleSet::validate] equalIntervalAngles3 attribute is not updated.";
    }

    if (oneSide_ && !distinguishOneSide()) {
        valid = false;
        lbWarn << "[SampleSet::validate] oneSide attribute is not updated.";
    }

    // Wavelengths
    if (wavelengths_.allFinite()) {
        if (wavelengths_.minCoeff() < 0.0f) {
            lbWarn
                << "[SampleSet::validate] The negative wavelength(s) is found:\n\t"
                << wavelengths_.format(LB_EIGEN_IO_FMT);
        }
        else {
            lbInfo << "[SampleSet::validate] Wavelengths are valid.";
        }
    }
    else {
        valid = false;
        lbWarn
            << "[SampleSet::validate] The invalid wavelength(s) is found:\n\t"
            << wavelengths_.format(LB_EIGEN_IO_FMT);
    }

    return valid;
}

void SampleSet::updateAngleAttributes()
{
    updateEqualIntervalAngles();
    updateOneSide();
}

void SampleSet::resizeAngles(int numAngles0,
                             int numAngles1,
                             int numAngles2,
                             int numAngles3)
{
    assert(numAngles0 > 0 && numAngles1 > 0 && numAngles2 > 0 && numAngles3 > 0);

    angles0_.resize(numAngles0);
    angles1_.resize(numAngles1);
    angles2_.resize(numAngles2);
    angles3_.resize(numAngles3);

    size_t numSamples = numAngles0 * numAngles1 * numAngles2 * numAngles3;
    spectra_.resize(numSamples);
}

void SampleSet::resizeWavelengths(int numWavelengths)
{
    assert(numWavelengths > 0);

    size_t numSamples = angles0_.size() * angles1_.size() * angles2_.size() * angles3_.size();

    for (size_t i = 0; i < numSamples; ++i) {
        spectra_.at(i) = Spectrum::Zero(numWavelengths);
    }

    wavelengths_.resize(numWavelengths);
}

void SampleSet::updateEqualIntervalAngles()
{
    equalIntervalAngles0_ = array_util::isEqualInterval(angles0_);
    equalIntervalAngles1_ = array_util::isEqualInterval(angles1_);
    equalIntervalAngles2_ = array_util::isEqualInterval(angles2_);
    equalIntervalAngles3_ = array_util::isEqualInterval(angles3_);

    lbInfo << "[SampleSet::updateEqualIntervalAngles] Angle0: " << equalIntervalAngles0_;
    lbInfo << "[SampleSet::updateEqualIntervalAngles] Angle1: " << equalIntervalAngles1_;
    lbInfo << "[SampleSet::updateEqualIntervalAngles] Angle2: " << equalIntervalAngles2_;
    lbInfo << "[SampleSet::updateEqualIntervalAngles] Angle3: " << equalIntervalAngles3_;
}

bool SampleSet::distinguishOneSide() const
{
    bool containing_0_PI = false;
    bool containing_PI_2PI = false;

    constexpr double offset = EPSILON_D * 2;

    for (int i = 0; i < angles3_.size(); ++i) {
        double angle = angles3_[i];

        if (angle > offset && angle < PI_D - offset * PI_D) {
            containing_0_PI = true;
        }

        if (angle > PI_D + offset * PI_D && angle < TAU_D - offset * TAU_D) {
            containing_PI_2PI = true;
        }
    }

    return (!containing_0_PI || !containing_PI_2PI);
}

void SampleSet::updateOneSide()
{
    oneSide_ = distinguishOneSide();

    lbInfo << "[SampleSet::updateOneSide] " << oneSide_;
}
