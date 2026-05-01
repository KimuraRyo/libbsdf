// =================================================================== //
// Copyright (C) 2026 Kimura Ryo                                       //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#include <libbsdf/Brdf/DistortedSphericalCoordinatesBrdf.h>

using namespace lb;

DistortedSphericalCoordinatesBrdf::DistortedSphericalCoordinatesBrdf(int        numInTheta,
                                                                     int        numInPhi,
                                                                     int        numDistTheta,
                                                                     int        numDistPhi,
                                                                     ColorModel colorModel,
                                                                     int        numWavelengths,
                                                                     bool       equalIntervalAngles)
    : BaseBrdf(numInTheta,
               numInPhi,
               numDistTheta,
               numDistPhi,
               colorModel,
               numWavelengths,
               equalIntervalAngles)
{
}

DistortedSphericalCoordinatesBrdf::DistortedSphericalCoordinatesBrdf(const Brdf&   brdf,
                                                                     const Arrayd& inThetaAngles,
                                                                     const Arrayd& inPhiAngles,
                                                                     const Arrayd& distThetaAngles,
                                                                     const Arrayd& distPhiAngles)
    : BaseBrdf(brdf, inThetaAngles, inPhiAngles, distThetaAngles, distPhiAngles)
{
}

DistortedSphericalCoordinatesBrdf::DistortedSphericalCoordinatesBrdf(const Brdf& brdf,
                                                                     int         numInTheta,
                                                                     int         numInPhi,
                                                                     int         numDistTheta,
                                                                     int         numDistPhi)
    : BaseBrdf(brdf, numInTheta, numInPhi, numDistTheta, numDistPhi)
{
}

DistortedSphericalCoordinatesBrdf::DistortedSphericalCoordinatesBrdf(int        numInTheta,
                                                                     int        numInPhi,
                                                                     int        numDistTheta,
                                                                     int        numDistPhi,
                                                                     double     distAngleExponent,
                                                                     ColorModel colorModel,
                                                                     int        numWavelengths,
                                                                     double     refractiveIndex)
    : BaseBrdf(numInTheta, numInPhi, numDistTheta, numDistPhi, colorModel, numWavelengths, true)
{
    Arrayd& distThetaAngles = samples_->getAngles2();
    distThetaAngles = array_util::createExponential<Arrayd>(
        static_cast<int>(distThetaAngles.size()), CoordSys::MAX_ANGLE2, distAngleExponent);

    if (numInPhi == 1) {
        // Adjust so that the density increases around an azimuth of 180 degrees.
        Arrayd& distPhiAngles = samples_->getAngles3();
        Arrayd  halfDistPhiAngles = array_util::createExponential<Arrayd>(
            static_cast<int>(distPhiAngles.size() / 2 + 1), PI_D, distAngleExponent);

        for (int i = 0; i < distPhiAngles.size() / 2; ++i) {
            distPhiAngles[i] = PI_D - halfDistPhiAngles[halfDistPhiAngles.size() - 1 - i];
            if (distPhiAngles.size() % 2 == 0) {
                distPhiAngles[distPhiAngles.size() / 2 + i] = PI_D + halfDistPhiAngles[i + 1];
            }
            else {
                distPhiAngles[distPhiAngles.size() / 2 + i] = PI_D + halfDistPhiAngles[i];
            }
        }
    }

    samples_->updateAngleAttributes();

    if (refractiveIndex != 1) {
        setupSpecularOffsets(refractiveIndex);
    }
}

DistortedSphericalCoordinatesBrdf::DistortedSphericalCoordinatesBrdf(
    const SphericalCoordinatesBrdf& brdf,
    int                             numDistTheta,
    int                             numDistPhi)
    : BaseBrdf()
{
    const SampleSet* ss = brdf.getSampleSet();

    Arrayd inThetaAngles = ss->getAngles0();
    Arrayd inPhiAngles = ss->getAngles1();

    if (inPhiAngles.size() == 1) {
        inPhiAngles[0] = 0;
    }

    Arrayd distThetaAngles =
        array_util::createExponential<Arrayd>(numDistTheta, CoordSys::MAX_ANGLE2, 2.0);

    Arrayd distPhiAngles;
    if (inPhiAngles.size() == 1) {
        // Adjust so that the density increases around an azimuth of 180 degrees.
        distPhiAngles.resize(numDistPhi);
        Arrayd halfDistPhiAngles = array_util::createExponential<Arrayd>(
            static_cast<int>(distPhiAngles.size() / 2 + 1), PI_D, 2.0);

        for (int i = 0; i < distPhiAngles.size() / 2; ++i) {
            distPhiAngles[i] = PI_D - halfDistPhiAngles[halfDistPhiAngles.size() - 1 - i];
            if (distPhiAngles.size() % 2 == 0) {
                distPhiAngles[distPhiAngles.size() / 2 + i] = PI_D + halfDistPhiAngles[i + 1];
            }
            else {
                distPhiAngles[distPhiAngles.size() / 2 + i] = PI_D + halfDistPhiAngles[i];
            }
        }
    }
    else {
        distPhiAngles = Arrayd::LinSpaced(numDistPhi, 0.0, CoordSys::MAX_ANGLE3);
    }

    samples_ = new SampleSet(
        static_cast<int>(inThetaAngles.size()), static_cast<int>(inPhiAngles.size()),
        static_cast<int>(distThetaAngles.size()), static_cast<int>(distPhiAngles.size()),
        ss->getColorModel(), ss->getNumWavelengths());
    samples_->getAngles0() = inThetaAngles;
    samples_->getAngles1() = inPhiAngles;
    samples_->getAngles2() = distThetaAngles;
    samples_->getAngles3() = distPhiAngles;
    samples_->getWavelengths() = ss->getWavelengths();

    samples_->updateAngleAttributes();
    initializeSpectra(brdf);

    reductionType_ = brdf.getReductionType();
    sourceType_ = brdf.getSourceType();
}

DistortedSphericalCoordinatesBrdf::DistortedSphericalCoordinatesBrdf(
    const DistortedSphericalCoordinatesBrdf& brdf)
    : BaseBrdf(brdf), specularOffsets_(brdf.getSpecularOffsets())
{
}

DistortedSphericalCoordinatesBrdf::~DistortedSphericalCoordinatesBrdf() {}

DistortedSphericalCoordinatesBrdf* DistortedSphericalCoordinatesBrdf::clone() const
{
    return new DistortedSphericalCoordinatesBrdf(*this);
}

void DistortedSphericalCoordinatesBrdf::setupSpecularOffsets(double refractiveIndex)
{
    for (int i = 0; i < getNumInTheta(); ++i) {
        double inTheta = getInTheta(i);
        double sinT = std::min(std::sin(inTheta) / refractiveIndex, 1.0);
        double refractedTheta = std::asin(sinT);
        double offset = refractedTheta - inTheta;

        setSpecularOffset(i, offset);
    }
}

bool DistortedSphericalCoordinatesBrdf::validate(bool verbose) const
{
    bool valid = BaseBrdf::validate(verbose);

    if (specularOffsets_.size() == 0)
        return valid;

    if (specularOffsets_.size() != getNumInTheta()) {
        valid = false;
        lbWarn
            << "[DistortedSphericalCoordinatesBrdf::validate] The number of specular offsets is invalid."
            << "\n\tSpecular offsets: " << specularOffsets_.size()
            << "\n\tIncoming polar angles: " << getNumInTheta();
    }
    else if (!specularOffsets_.allFinite()) {
        valid = false;
        lbWarn << "[DistortedSphericalCoordinatesBrdf::validate] The array of specular offset is invalid.";
    }
    else if (specularOffsets_.minCoeff() < -CoordSys::MAX_ANGLE0 ||
             specularOffsets_.maxCoeff() >  CoordSys::MAX_ANGLE0) {
        valid = false;
        lbWarn << "[DistortedSphericalCoordinatesBrdf::validate] The angle(s) in specular offsets is outside of range.";
    }
    else {
        lbInfo << "[DistortedSphericalCoordinatesBrdf::validate] The array of specular offset is valid.";
    }

    return valid;
}

bool DistortedSphericalCoordinatesBrdf::expandAngles(bool angle0Expanded,
                                                     bool angle1Expanded,
                                                     bool angle2Expanded,
                                                     bool angle3Expanded)
{
    const Arrayd& angles0 = samples_->getAngles0();
    Arrayd        origAngles0 = angles0;

    bool expanded =
        BaseBrdf::expandAngles(angle0Expanded, angle1Expanded, angle2Expanded, angle3Expanded);

    if (specularOffsets_.size() == 0)
        return expanded;

    if (origAngles0.size() == angles0.size())
        return expanded;

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