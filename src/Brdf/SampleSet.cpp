// =================================================================== //
// Copyright (C) 2014-2015 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#include <libbsdf/Brdf/SampleSet.h>

#include <algorithm>
#include <iostream>

#include <libbsdf/Common/Utility.h>
#include <libbsdf/Common/SpectrumUtility.h>

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
                       equalIntervalAngles3_(false)
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

void SampleSet::zeroOutSamples()
{
    angles0_ = Arrayf::Zero(angles0_.size());
    angles1_ = Arrayf::Zero(angles1_.size());
    angles2_ = Arrayf::Zero(angles2_.size());
    angles3_ = Arrayf::Zero(angles3_.size());

    for (auto it = spectra_.begin(); it != spectra_.end(); ++it) {
        *it = Spectrum::Zero(it->size());
    }

    wavelengths_ = Arrayf::Zero(wavelengths_.size());
}

void SampleSet::checkEqualIntervalAngles()
{
    equalIntervalAngles0_ = isEqualInterval(angles0_);
    equalIntervalAngles1_ = isEqualInterval(angles1_);
    equalIntervalAngles2_ = isEqualInterval(angles2_);
    equalIntervalAngles3_ = isEqualInterval(angles3_);

    std::cout << "[SampleSet::checkEqualIntervalAngles] Angle0: " << equalIntervalAngles0_ << std::endl;
    std::cout << "[SampleSet::checkEqualIntervalAngles] Angle1: " << equalIntervalAngles1_ << std::endl;
    std::cout << "[SampleSet::checkEqualIntervalAngles] Angle2: " << equalIntervalAngles2_ << std::endl;
    std::cout << "[SampleSet::checkEqualIntervalAngles] Angle3: " << equalIntervalAngles3_ << std::endl;
}

void SampleSet::resizeAngles(int numAngles0,
                             int numAngles1,
                             int numAngles2,
                             int numAngles3)
{
    assert(numAngles0 > 0 && numAngles1 > 0 && numAngles2 > 0 && numAngles3 > 0);

    numAngles0_ = numAngles0;
    numAngles1_ = numAngles1;
    numAngles2_ = numAngles2;
    numAngles3_ = numAngles3;

    int numSamples = numAngles0 * numAngles1 * numAngles2 * numAngles3;
    spectra_.resize(numSamples);

    angles0_.resize(numAngles0);
    angles1_.resize(numAngles1);
    angles2_.resize(numAngles2);
    angles3_.resize(numAngles3);
}

void SampleSet::resizeWavelengths(int numWavelengths)
{
    assert(numWavelengths > 0);

    int numSamples = numAngles0_ * numAngles1_ * numAngles2_ * numAngles3_;

    for (int i = 0; i < numSamples; ++i) {
        Spectrum spectrum;
        spectrum.resize(numWavelengths);
        spectra_.at(i) = spectrum;
    }

    wavelengths_.resize(numWavelengths);
}
