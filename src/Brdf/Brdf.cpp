// =================================================================== //
// Copyright (C) 2014-2019 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#include <libbsdf/Brdf/Brdf.h>

using namespace lb;

Brdf::Brdf(int          numAngles0,
           int          numAngles1,
           int          numAngles2,
           int          numAngles3,
           ColorModel   colorModel,
           int          numWavelengths)
           : samples_(new SampleSet(numAngles0,
                                    numAngles1,
                                    numAngles2,
                                    numAngles3,
                                    colorModel,
                                    numWavelengths)),
             sourceType_(UNKNOWN_SOURCE)
{
    lbTrace << "[Brdf::Brdf]";
}

Brdf::Brdf() : samples_(0),
               sourceType_(UNKNOWN_SOURCE)
{
    lbTrace << "[Brdf::Brdf]";
}

Brdf::Brdf(const Brdf& brdf) : samples_(new SampleSet(*brdf.getSampleSet())),
                               sourceType_(brdf.getSourceType())
{
    lbTrace << "[Brdf::Brdf]";
}

Brdf::~Brdf()
{
    lbTrace << "[Brdf::~Brdf]";
    delete samples_;
}
