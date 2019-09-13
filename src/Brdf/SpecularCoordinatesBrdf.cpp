// =================================================================== //
// Copyright (C) 2014-2019 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#include <libbsdf/Brdf/SpecularCoordinatesBrdf.h>

#include <libbsdf/Common/Log.h>

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
    Arrayf& specThetaAngles = samples_->getAngles2();
    specThetaAngles = createExponentialArray<Arrayf>(static_cast<int>(specThetaAngles.size()),
                                                     CoordSys::MAX_ANGLE2,
                                                     specThetaExponent);

    samples_->updateAngleAttributes();

    if (refractiveIndex != 1.0f) {
        for (int i = 0; i < getNumInTheta(); ++i) {
            float inTheta = getInTheta(i);
            float sinT = std::min(std::sin(inTheta) / refractiveIndex, 1.0f);
            float refractedTheta = std::asin(sinT);
            setSpecularOffset(i, refractedTheta - inTheta);
        }
    }
}

SpecularCoordinatesBrdf::SpecularCoordinatesBrdf(const SphericalCoordinatesBrdf& brdf,
                                                 int                             numSpecTheta,
                                                 int                             numSpecPhi)
                                                 : BaseBrdf()
{
    const SampleSet* ss = brdf.getSampleSet();

    Arrayf inThetaAngles   = ss->getAngles0();
    Arrayf inPhiAngles     = ss->getAngles1();

    if (inPhiAngles.size() == 1) {
        inPhiAngles[0] = 0.0f;
    }

    Arrayf specThetaAngles = createExponentialArray<Arrayf>(numSpecTheta, CoordSys::MAX_ANGLE2, 2.0f);
    Arrayf specPhiAngles   = Arrayf::LinSpaced(numSpecPhi, 0.0, CoordSys::MAX_ANGLE3);

    samples_ = new SampleSet(static_cast<int>(inThetaAngles.size()),
                             static_cast<int>(inPhiAngles.size()),
                             static_cast<int>(specThetaAngles.size()),
                             static_cast<int>(specPhiAngles.size()),
                             ss->getColorModel(),
                             ss->getNumWavelengths());
    samples_->getAngles0() = inThetaAngles;
    samples_->getAngles1() = inPhiAngles;
    samples_->getAngles2() = specThetaAngles;
    samples_->getAngles3() = specPhiAngles;
    samples_->getWavelengths() = ss->getWavelengths();

    samples_->updateAngleAttributes();
    initializeSpectra(brdf);

    sourceType_ = brdf.getSourceType();
}

SpecularCoordinatesBrdf::SpecularCoordinatesBrdf(const SpecularCoordinatesBrdf& brdf)
                                                 : BaseBrdf(brdf),
                                                   specularOffsets_(brdf.getSpecularOffsets()) {}

SpecularCoordinatesBrdf::~SpecularCoordinatesBrdf() {}

SpecularCoordinatesBrdf* SpecularCoordinatesBrdf::clone() const
{
    return new SpecularCoordinatesBrdf(*this);
}

bool SpecularCoordinatesBrdf::validate(bool verbose) const
{
    bool valid = BaseBrdf::validate(verbose);

    if (specularOffsets_.size() == 0) return valid;

    if (specularOffsets_.size() != getNumInTheta()) {
        valid = false;
        lbWarn
            << "[SpecularCoordinatesBrdf::validate] The number of specular offsets is invalid."
            << "\n\tSpecular offsets: " << specularOffsets_.size()
            << "\n\tIncoming polar angles: " << getNumInTheta();
    }
    else if (!specularOffsets_.allFinite()) {
        valid = false;
        lbWarn << "[SpecularCoordinatesBrdf::validate] The array of specular offset is invalid.";
    }
    else if (specularOffsets_.minCoeff() < -CoordSys::MAX_ANGLE0 ||
             specularOffsets_.maxCoeff() >  CoordSys::MAX_ANGLE0) {
        valid = false;
        lbWarn << "[SpecularCoordinatesBrdf::validate] The angle(s) in specular offsets is outside of range.";
    }
    else {
        lbInfo << "[SpecularCoordinatesBrdf::validate] The array of specular offset is valid.";
    }

    return valid;
}
