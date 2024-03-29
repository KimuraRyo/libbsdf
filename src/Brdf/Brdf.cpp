// =================================================================== //
// Copyright (C) 2014-2020 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#include <libbsdf/Brdf/Brdf.h>

#include <libbsdf/Brdf/Initializer.h>
#include <libbsdf/Brdf/LinearInterpolator.h>

using namespace lb;

Brdf::Brdf(int        numAngles0,
           int        numAngles1,
           int        numAngles2,
           int        numAngles3,
           ColorModel colorModel,
           int        numWavelengths)
    : samples_(new SampleSet(numAngles0,
                             numAngles1,
                             numAngles2,
                             numAngles3,
                             colorModel,
                             numWavelengths)),
      reductionType_(ReductionType::NONE),
      sourceType_(UNKNOWN_SOURCE)
{
    lbTrace << "[Brdf::Brdf]";
}

Brdf::Brdf() : samples_(0), reductionType_(ReductionType::NONE), sourceType_(UNKNOWN_SOURCE)
{
    lbTrace << "[Brdf::Brdf]";
}

Brdf::Brdf(const Brdf& brdf)
    : samples_(new SampleSet(*brdf.getSampleSet())),
      reductionType_(brdf.getReductionType()),
      sourceType_(brdf.getSourceType())
{
    lbTrace << "[Brdf::Brdf]";
}

Brdf::~Brdf()
{
    lbTrace << "[Brdf::~Brdf] " << name_;
    delete samples_;
}

void Brdf::initializeSpectra(const Brdf& brdf)
{
    lb::initializeSpectra<LinearInterpolator>(brdf, this);
}

void Brdf::setName(const std::string& name)
{
    lbTrace << "[Brdf::setName] " << name;
    name_ = name;
}
