// =================================================================== //
// Copyright (C) 2014-2016 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#include <libbsdf/Brdf/HalfDifferenceCoordinatesBrdf.h>

using namespace lb;

HalfDifferenceCoordinatesBrdf::HalfDifferenceCoordinatesBrdf(int        numHalfTheta,
                                                             int        numHalfPhi,
                                                             int        numDiffTheta,
                                                             int        numDiffPhi,
                                                             ColorModel colorModel,
                                                             int        numWavelengths,
                                                             bool       equalIntervalAngles)
                                                             : BaseBrdf(numHalfTheta,
                                                                        numHalfPhi,
                                                                        numDiffTheta,
                                                                        numDiffPhi,
                                                                        colorModel,
                                                                        numWavelengths,
                                                                        equalIntervalAngles) {}

HalfDifferenceCoordinatesBrdf::HalfDifferenceCoordinatesBrdf(const Brdf&    brdf,
                                                             const Arrayf&  halfThetaAngles,
                                                             const Arrayf&  halfPhiAngles,
                                                             const Arrayf&  diffThetaAngles,
                                                             const Arrayf&  diffPhiAngles)
                                                             : BaseBrdf(brdf,
                                                                        halfThetaAngles,
                                                                        halfPhiAngles,
                                                                        diffThetaAngles,
                                                                        diffPhiAngles) {}

HalfDifferenceCoordinatesBrdf::HalfDifferenceCoordinatesBrdf(const Brdf&    brdf,
                                                             int            numHalfTheta,
                                                             int            numHalfPhi,
                                                             int            numDiffTheta,
                                                             int            numDiffPhi)
                                                             : BaseBrdf(brdf,
                                                                        numHalfTheta,
                                                                        numHalfPhi,
                                                                        numDiffTheta,
                                                                        numDiffPhi) {}

HalfDifferenceCoordinatesBrdf::HalfDifferenceCoordinatesBrdf(const HalfDifferenceCoordinatesBrdf& brdf)
                                                             : BaseBrdf(brdf) {}

HalfDifferenceCoordinatesBrdf::~HalfDifferenceCoordinatesBrdf() {}

HalfDifferenceCoordinatesBrdf* HalfDifferenceCoordinatesBrdf::clone() const
{
    return new HalfDifferenceCoordinatesBrdf(*this);
}
