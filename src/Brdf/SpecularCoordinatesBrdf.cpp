// =================================================================== //
// Copyright (C) 2014-2019 Kimura Ryo                                  //
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

SpecularCoordinatesBrdf::SpecularCoordinatesBrdf(int        numInTheta,
                                                 int        numInPhi,
                                                 int        numSpecTheta,
                                                 int        numSpecPhi,
                                                 float      specThetaExponent,
                                                 ColorModel colorModel,
                                                 int        numWavelengths,
                                                 float      refractiveIndex)
                                                 : BaseBrdf(numInTheta,
                                                            numInPhi,
                                                            numSpecTheta,
                                                            numSpecPhi,
                                                            colorModel,
                                                            numWavelengths,
                                                            true)
{
    // Create narrow intervals near specular directions.
    Arrayf& specThetaAngles = samples_->getAngles2();
    for (int i = 1; i < specThetaAngles.size() - 1; ++i) {
        Arrayf::Scalar ratio = specThetaAngles[i] / SpecularCoordinateSystem::MAX_ANGLE2;
        ratio = std::pow(ratio, static_cast<Arrayf::Scalar>(specThetaExponent));
        specThetaAngles[i] = ratio * SpecularCoordinateSystem::MAX_ANGLE2;
    }

    if (refractiveIndex != 1.0f) {
        for (int i = 0; i < getNumInTheta(); ++i) {
            float inTheta = getInTheta(i);
            float sinT = std::min(std::sin(inTheta) / refractiveIndex, 1.0f);
            float refractedTheta = std::asin(sinT);
            setSpecularOffset(i, refractedTheta - inTheta);
        }
    }
}

SpecularCoordinatesBrdf::SpecularCoordinatesBrdf(const SpecularCoordinatesBrdf& brdf)
                                                 : BaseBrdf(brdf),
                                                   specularOffsets_(brdf.getSpecularOffsets()) {}

SpecularCoordinatesBrdf::~SpecularCoordinatesBrdf() {}

SpecularCoordinatesBrdf* SpecularCoordinatesBrdf::clone() const
{
    return new SpecularCoordinatesBrdf(*this);
}
