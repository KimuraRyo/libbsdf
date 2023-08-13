// =================================================================== //
// Copyright (C) 2014-2023 Kimura Ryo                                  //
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
               equalIntervalAngles)
{
}

SpecularCoordinatesBrdf::SpecularCoordinatesBrdf(const Brdf&   brdf,
                                                 const Arrayd& inThetaAngles,
                                                 const Arrayd& inPhiAngles,
                                                 const Arrayd& specThetaAngles,
                                                 const Arrayd& specPhiAngles)
    : BaseBrdf(brdf, inThetaAngles, inPhiAngles, specThetaAngles, specPhiAngles)
{
}

SpecularCoordinatesBrdf::SpecularCoordinatesBrdf(const Brdf& brdf,
                                                 int         numInTheta,
                                                 int         numInPhi,
                                                 int         numSpecTheta,
                                                 int         numSpecPhi)
    : BaseBrdf(brdf, numInTheta, numInPhi, numSpecTheta, numSpecPhi)
{
}

SpecularCoordinatesBrdf::SpecularCoordinatesBrdf(int        numInTheta,
                                                 int        numInPhi,
                                                 int        numSpecTheta,
                                                 int        numSpecPhi,
                                                 double     specThetaExponent,
                                                 ColorModel colorModel,
                                                 int        numWavelengths,
                                                 double     refractiveIndex)
    : BaseBrdf(numInTheta, numInPhi, numSpecTheta, numSpecPhi, colorModel, numWavelengths, true)
{
    Arrayd& specThetaAngles = samples_->getAngles2();
    specThetaAngles = array_util::createExponential<Arrayd>(
        static_cast<int>(specThetaAngles.size()), CoordSys::MAX_ANGLE2, specThetaExponent);

    samples_->updateAngleAttributes();

    if (refractiveIndex != 1) {
        setupSpecularOffsets(refractiveIndex);
    }
}

SpecularCoordinatesBrdf::SpecularCoordinatesBrdf(const SphericalCoordinatesBrdf& brdf,
                                                 int                             numSpecTheta,
                                                 int                             numSpecPhi)
    : BaseBrdf()
{
    const SampleSet* ss = brdf.getSampleSet();

    Arrayd inThetaAngles = ss->getAngles0();
    Arrayd inPhiAngles = ss->getAngles1();

    if (inPhiAngles.size() == 1) {
        inPhiAngles[0] = 0;
    }

    Arrayd specThetaAngles =
        array_util::createExponential<Arrayd>(numSpecTheta, CoordSys::MAX_ANGLE2, 2.0);
    Arrayd specPhiAngles = Arrayd::LinSpaced(numSpecPhi, 0.0, CoordSys::MAX_ANGLE3);

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

    reductionType_ = brdf.getReductionType();
    sourceType_ = brdf.getSourceType();
}

SpecularCoordinatesBrdf::SpecularCoordinatesBrdf(const SpecularCoordinatesBrdf& brdf)
    : BaseBrdf(brdf), specularOffsets_(brdf.getSpecularOffsets())
{
}

SpecularCoordinatesBrdf::~SpecularCoordinatesBrdf() {}

SpecularCoordinatesBrdf* SpecularCoordinatesBrdf::clone() const
{
    return new SpecularCoordinatesBrdf(*this);
}

void SpecularCoordinatesBrdf::setupSpecularOffsets(double refractiveIndex)
{
    for (int i = 0; i < getNumInTheta(); ++i) {
        double inTheta = getInTheta(i);
        double sinT = std::min(std::sin(inTheta) / refractiveIndex, 1.0);
        double refractedTheta = std::asin(sinT);
        double offset = refractedTheta - inTheta;

        setSpecularOffset(i, offset);
    }
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

bool SpecularCoordinatesBrdf::expandAngles(bool angle0Expanded,
                                           bool angle1Expanded,
                                           bool angle2Expanded,
                                           bool angle3Expanded)
{
    const Arrayd& angles0 = samples_->getAngles0();
    Arrayd        origAngles0 = angles0;

    bool expanded = BaseBrdf::expandAngles(angle0Expanded,
                                           angle1Expanded,
                                           angle2Expanded,
                                           angle3Expanded);

    if (specularOffsets_.size() == 0) return expanded;

    if (origAngles0.size() == angles0.size()) return expanded;

    // Add the first element.
    if (origAngles0[0] != angles0[0]) {
        Arrayd origOffsets = specularOffsets_;

        specularOffsets_.resize(origOffsets.size() + 1);
        specularOffsets_[0] = origOffsets[0];
        for (int i = 1; i < specularOffsets_.size(); ++i) {
            specularOffsets_[i] = origOffsets[i - 1];
        }
    }

    // Add the last element.
    if (origAngles0[origAngles0.size() - 1] != angles0[angles0.size() - 1]) {
        Arrayd origOffsets = specularOffsets_;

        specularOffsets_.resize(origOffsets.size() + 1);
        for (int i = 0; i < specularOffsets_.size() - 1; ++i) {
            specularOffsets_[i] = origOffsets[i];
        }
        specularOffsets_[specularOffsets_.size() - 1] = origOffsets[origOffsets.size() - 1];
    }

    return expanded;
}
