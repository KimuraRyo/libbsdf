// =================================================================== //
// Copyright (C) 2014-2023 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#include <libbsdf/Brdf/SphericalCoordinatesBrdf.h>

using namespace lb;

SphericalCoordinatesBrdf::SphericalCoordinatesBrdf(int        numInTheta,
                                                   int        numInPhi,
                                                   int        numOutTheta,
                                                   int        numOutPhi,
                                                   ColorModel colorModel,
                                                   int        numWavelengths,
                                                   bool       equalIntervalAngles)
    : BaseBrdf(numInTheta,
               numInPhi,
               numOutTheta,
               numOutPhi,
               colorModel,
               numWavelengths,
               equalIntervalAngles)
{
}

SphericalCoordinatesBrdf::SphericalCoordinatesBrdf(const Brdf&   brdf,
                                                   const Arrayd& inThetaAngles,
                                                   const Arrayd& inPhiAngles,
                                                   const Arrayd& outThetaAngles,
                                                   const Arrayd& outPhiAngles)
    : BaseBrdf(brdf, inThetaAngles, inPhiAngles, outThetaAngles, outPhiAngles)
{
}

SphericalCoordinatesBrdf::SphericalCoordinatesBrdf(const Brdf& brdf,
                                                   int         numInTheta,
                                                   int         numInPhi,
                                                   int         numOutTheta,
                                                   int         numOutPhi)
    : BaseBrdf(brdf, numInTheta, numInPhi, numOutTheta, numOutPhi)
{
}

SphericalCoordinatesBrdf::SphericalCoordinatesBrdf(const SphericalCoordinatesBrdf& brdf)
    : BaseBrdf(brdf)
{
}

SphericalCoordinatesBrdf::~SphericalCoordinatesBrdf() {}

SphericalCoordinatesBrdf* SphericalCoordinatesBrdf::clone() const
{
    return new SphericalCoordinatesBrdf(*this);
}
