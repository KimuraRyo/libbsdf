// =================================================================== //
// Copyright (C) 2014-2015 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#include <libbsdf/Brdf/SpecularCoordinatesBrdf.h>

using namespace lb;

SpecularCoordinatesBrdf::SpecularCoordinatesBrdf(int        numInTheta,
                                                 int        numInPhi,
                                                 int        numSpecTheta,
                                                 int        numSpecPhi,
                                                 ColorModel colorModel,
                                                 int        numWavelengths,
                                                 bool       equalIntervalAngles)
                                                 : BaseBrdf(numInTheta,
                                                            numInPhi,
                                                            numSpecTheta,
                                                            numSpecPhi,
                                                            colorModel,
                                                            numWavelengths,
                                                            equalIntervalAngles) {}

SpecularCoordinatesBrdf::SpecularCoordinatesBrdf(const Brdf&    brdf,
                                                 const Arrayf&  inThetaAngles,
                                                 const Arrayf&  inPhiAngles,
                                                 const Arrayf&  specThetaAngles,
                                                 const Arrayf&  specPhiAngles)
                                                 : BaseBrdf(brdf,
                                                            inThetaAngles,
                                                            inPhiAngles,
                                                            specThetaAngles,
                                                            specPhiAngles) {}

SpecularCoordinatesBrdf::SpecularCoordinatesBrdf(const Brdf&    brdf,
                                                 int            numInTheta,
                                                 int            numInPhi,
                                                 int            numSpecTheta,
                                                 int            numSpecPhi)
                                                 : BaseBrdf(brdf,
                                                            numInTheta,
                                                            numInPhi,
                                                            numSpecTheta,
                                                            numSpecPhi) {}

SpecularCoordinatesBrdf::SpecularCoordinatesBrdf(const SpecularCoordinatesBrdf& brdf)
                                                 : BaseBrdf(brdf) {}

SpecularCoordinatesBrdf::~SpecularCoordinatesBrdf() {}
