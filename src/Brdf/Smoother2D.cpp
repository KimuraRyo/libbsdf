// =================================================================== //
// Copyright (C) 2017-2019 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#include <libbsdf/Brdf/Smoother2D.h>

#include <libbsdf/Brdf/CatmullRomSplineInterpolator.h>
#include <libbsdf/Brdf/Initializer.h>
#include <libbsdf/Brdf/LinearInterpolator.h>
#include <libbsdf/Brdf/Sampler.h>

#include <libbsdf/Common/Array.h>

using namespace lb;

Smoother2D::Smoother2D(SampleSet2D* samples)
                       : samples_(samples),
                         diffThreshold_(0.001f),
                         maxIteration0_(2),
                         maxIteration1_(2),
                         minAngleInterval_(toRadian(0.1f)) {}

void Smoother2D::smooth()
{
    initializeAngles();

    for (int i = 0; i < maxIteration0_; ++i) {
        if (!insertAngle0()) {
            break;
        }
        updateSamples();
    }

    for (int i = 0; i < maxIteration1_; ++i) {
        if (!insertAngle1()) {
            break;
        }
        updateSamples();
    }
}

void Smoother2D::initializeAngles()
{
    for (int i = 0; i < samples_->getNumTheta(); ++i) {
        angles0_.insert(samples_->getTheta(i));
    }

    for (int i = 0; i < samples_->getNumPhi(); ++i) {
        angles1_.insert(samples_->getPhi(i));
    }
}

bool Smoother2D::insertAngle0()
{
    if (samples_->getNumTheta() < 4) {
        return false;
    }

    bool inserted = false;

    for (int i1 = 0; i1 < samples_->getNumPhi(); ++i1) {
        for (int i0 = 1; i0 < samples_->getNumTheta() - 2; ++i0) {
            Vec2f angles(samples_->getTheta(i0),
                         samples_->getPhi(i1));

            Vec2f nextAngles(samples_->getTheta(i0 + 1),
                             samples_->getPhi(i1));

            if (insertAngle(angles0_, 0, angles, nextAngles)) {
                inserted = true;
            }
        }
    }

    return inserted;
}

bool Smoother2D::insertAngle1()
{
    if (samples_->getNumPhi() < 4) {
        return false;
    }

    bool inserted = false;

    for (int i0 = 0; i0 < samples_->getNumTheta(); ++i0) {
        for (int i1 = 0; i1 < samples_->getNumPhi() - 1; ++i1) {
            Vec2f angles(samples_->getTheta(i0),
                         samples_->getPhi(i1));

            Vec2f nextAngles(samples_->getTheta(i0),
                             samples_->getPhi(i1 + 1));

            if (insertAngle(angles1_, 1, angles, nextAngles)) {
                inserted = true;
            }
        }
    }

    return inserted;
}

bool Smoother2D::insertAngle(std::set<Arrayf::Scalar>&  angleSet,
                             int                        angleSuffix,
                             const Vec2f&               angles,
                             const Vec2f&               nextAngles)
{
    if (std::abs(angles[angleSuffix] - nextAngles[angleSuffix]) < minAngleInterval_) {
        return false;
    }

    Vec2f midAngles = (angles + nextAngles) / 2.0f;

    Spectrum lerpSp;
    LinearInterpolator::getSpectrum(*samples_, midAngles[0], midAngles[1], &lerpSp);

    Spectrum crSp;
    CatmullRomSplineInterpolator::getSpectrum(*samples_, midAngles[0], midAngles[1], &crSp);

    Spectrum diffSp = (lerpSp - crSp).abs();
    if (diffSp.maxCoeff() > diffThreshold_) {
        auto ret = angleSet.insert(midAngles[angleSuffix]);
        if (ret.second) {
            return true;
        }
    }

    return false;
}

void Smoother2D::updateSamples()
{
    SampleSet2D* origSs2 = new SampleSet2D(*samples_);
    samples_->resizeAngles(static_cast<int>(angles0_.size()), static_cast<int>(angles1_.size()));

    copyArray(angles0_, &samples_->getThetaArray());
    copyArray(angles1_, &samples_->getPhiArray());
    samples_->updateAngleAttributes();

    lb::initializeSpectra<CatmullRomSplineInterpolator>(*origSs2, samples_);

    delete origSs2;
}
