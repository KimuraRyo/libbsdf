// =================================================================== //
// Copyright (C) 2017-2020 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#include <libbsdf/Brdf/Smoother.h>

#include <libbsdf/Brdf/CatmullRomSplineInterpolator.h>
#include <libbsdf/Brdf/Initializer.h>
#include <libbsdf/Brdf/LinearInterpolator.h>
#include <libbsdf/Brdf/Sampler.h>
#include <libbsdf/Brdf/SpecularCoordinatesBrdf.h>

#include <libbsdf/Common/Array.h>

using namespace lb;

Smoother::Smoother(Brdf* brdf)
                   : brdf_(brdf),
                     diffThreshold_(0.001f),
                     maxIteration0_(2),
                     maxIteration1_(2),
                     maxIteration2_(2),
                     maxIteration3_(2),
                     minAngleInterval_(toRadian(0.1f)),
                     specularPolarRegion_(0.0f) {}

void Smoother::smooth()
{
    initializeAngles();

    for (int i = 0; i < maxIteration0_; ++i) {
        if (!insertAngle0()) {
            break;
        }
        updateBrdf();
    }

    for (int i = 0; i < maxIteration1_; ++i) {
        if (!insertAngle1()) {
            break;
        }
        updateBrdf();
    }

    for (int i = 0; i < maxIteration2_; ++i) {
        if (!insertAngle2()) {
            break;
        }
        updateBrdf();
    }

    for (int i = 0; i < maxIteration3_; ++i) {
        if (!insertAngle3()) {
            break;
        }
        updateBrdf();
    }
}

void Smoother::initializeAngles()
{
    const SampleSet* ss = brdf_->getSampleSet();

    for (int i = 0; i < ss->getNumAngles0(); ++i) {
        angles0_.insert(ss->getAngle0(i));
    }

    for (int i = 0; i < ss->getNumAngles1(); ++i) {
        angles1_.insert(ss->getAngle1(i));
    }

    for (int i = 0; i < ss->getNumAngles2(); ++i) {
        angles2_.insert(ss->getAngle2(i));
    }

    for (int i = 0; i < ss->getNumAngles3(); ++i) {
        angles3_.insert(ss->getAngle3(i));
    }
}

bool Smoother::insertAngle0()
{
    const SampleSet* ss = brdf_->getSampleSet();

    if (ss->getNumAngles0() < 4) {
        return false;
    }

    bool inserted = false;

    for (int i1 = 0; i1 < ss->getNumAngles1(); ++i1) {
    for (int i2 = 0; i2 < ss->getNumAngles2(); ++i2) {
    for (int i3 = 0; i3 < ss->getNumAngles3(); ++i3) {
        for (int i0 = 1; i0 < ss->getNumAngles0() - 2; ++i0) {
            Vec4f angles(ss->getAngle0(i0),
                         ss->getAngle1(i1),
                         ss->getAngle2(i2),
                         ss->getAngle3(i3));

            Vec4f nextAngles(ss->getAngle0(i0 + 1),
                             ss->getAngle1(i1),
                             ss->getAngle2(i2),
                             ss->getAngle3(i3));

            if (insertAngle(angles0_, 0, angles, nextAngles)) {
                inserted = true;
            }
        }
    }}}

    return inserted;
}

bool Smoother::insertAngle1()
{
    const SampleSet* ss = brdf_->getSampleSet();

    if (ss->getNumAngles1() < 4) {
        return false;
    }

    bool inserted = false;

    for (int i0 = 0; i0 < ss->getNumAngles0(); ++i0) {
    for (int i2 = 0; i2 < ss->getNumAngles2(); ++i2) {
    for (int i3 = 0; i3 < ss->getNumAngles3(); ++i3) {
        for (int i1 = 0; i1 < ss->getNumAngles1() - 1; ++i1) {
            Vec4f angles(ss->getAngle0(i0),
                         ss->getAngle1(i1),
                         ss->getAngle2(i2),
                         ss->getAngle3(i3));

            Vec4f nextAngles(ss->getAngle0(i0),
                             ss->getAngle1(i1 + 1),
                             ss->getAngle2(i2),
                             ss->getAngle3(i3));

            if (insertAngle(angles1_, 1, angles, nextAngles)) {
                inserted = true;
            }
        }
    }}}

    return inserted;
}

bool Smoother::insertAngle2()
{
    const SampleSet* ss = brdf_->getSampleSet();

    if (ss->getNumAngles2() < 4) {
        return false;
    }

    bool inserted = false;

    for (int i0 = 0; i0 < ss->getNumAngles0(); ++i0) {
    for (int i1 = 0; i1 < ss->getNumAngles1(); ++i1) {
    for (int i3 = 0; i3 < ss->getNumAngles3(); ++i3) {
        for (int i2 = 1; i2 < ss->getNumAngles2() - 2; ++i2) {
            Vec4f angles(ss->getAngle0(i0),
                         ss->getAngle1(i1),
                         ss->getAngle2(i2),
                         ss->getAngle3(i3));

            Vec4f nextAngles(ss->getAngle0(i0),
                             ss->getAngle1(i1),
                             ss->getAngle2(i2 + 1),
                             ss->getAngle3(i3));

            if (dynamic_cast<SpecularCoordinatesBrdf*>(brdf_) &&
                angles[2] <= specularPolarRegion_) {
                continue;
            }


            if (insertAngle(angles2_, 2, angles, nextAngles)) {
                inserted = true;
            }
        }
    }}}

    return inserted;
}

bool Smoother::insertAngle3()
{
    const SampleSet* ss = brdf_->getSampleSet();

    if (ss->getNumAngles3() < 4) {
        return false;
    }

    bool inserted = false;

    for (int i0 = 0; i0 < ss->getNumAngles0(); ++i0) {
    for (int i1 = 0; i1 < ss->getNumAngles1(); ++i1) {
    for (int i2 = 0; i2 < ss->getNumAngles2(); ++i2) {
        for (int i3 = 0; i3 < ss->getNumAngles3() - 1; ++i3) {
            Vec4f angles(ss->getAngle0(i0),
                         ss->getAngle1(i1),
                         ss->getAngle2(i2),
                         ss->getAngle3(i3));

            Vec4f nextAngles(ss->getAngle0(i0),
                             ss->getAngle1(i1),
                             ss->getAngle2(i2),
                             ss->getAngle3(i3 + 1));

            if (insertAngle(angles3_, 3, angles, nextAngles)) {
                inserted = true;
            }
        }
    }}}

    return inserted;
}

bool Smoother::insertAngle(std::set<Arrayf::Scalar>&    angleSet,
                           int                          angleSuffix,
                           const Vec4f&                 angles,
                           const Vec4f&                 nextAngles)
{
    if (std::abs(angles[angleSuffix] - nextAngles[angleSuffix]) < minAngleInterval_) {
        return false;
    }

    Vec4f midAngles = (angles + nextAngles) / 2.0f;

    lb::Vec3 inDir, outDir;
    brdf_->toXyz(midAngles[0], midAngles[1], midAngles[2], midAngles[3], &inDir, &outDir);

    if (inDir.z() <= 0.0 || outDir.z() <= 0.0) {
        return false;
    }

    Spectrum lerpSp = Sampler::getSpectrum<LinearInterpolator>(*brdf_, inDir, outDir);
    Spectrum crSp   = Sampler::getSpectrum<CatmullRomSplineInterpolator>(*brdf_, inDir, outDir);

    Spectrum diffSp = (lerpSp - crSp).abs();
    if (diffSp.maxCoeff() > diffThreshold_) {
        auto ret = angleSet.insert(midAngles[angleSuffix]);
        if (ret.second) {
            return true;
        }
    }

    return false;
}

void Smoother::updateBrdf()
{
    Brdf* origBrdf = brdf_->clone();

    SampleSet* ss = brdf_->getSampleSet();
    ss->resizeAngles(static_cast<int>(angles0_.size()),
                     static_cast<int>(angles1_.size()),
                     static_cast<int>(angles2_.size()),
                     static_cast<int>(angles3_.size()));

    array_util::copy(angles0_, &ss->getAngles0());
    array_util::copy(angles1_, &ss->getAngles1());
    array_util::copy(angles2_, &ss->getAngles2());
    array_util::copy(angles3_, &ss->getAngles3());
    ss->updateAngleAttributes();

    initializeSpectra<CatmullRomSplineInterpolator>(*origBrdf, brdf_);

    delete origBrdf;
}
