// =================================================================== //
// Copyright (C) 2020-2023 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#include <libbsdf/Brdf/Optimizer.h>

#include <libbsdf/Brdf/Initializer.h>
#include <libbsdf/Brdf/LinearInterpolator.h>

#include <libbsdf/Common/Array.h>

using namespace lb;

Optimizer::Optimizer(Brdf* brdf, double diffThreshold, double ratioThreshold)
    : brdf_(brdf), diffThreshold_(diffThreshold), ratioThreshold_(ratioThreshold)
{
}

void Optimizer::optimize()
{
    setupAngles0();
    setupAngles1();
    setupAngles2();
    setupAngles3();

    updateBrdf();
}

void Optimizer::setupAngles0()
{
    const SampleSet* ss = brdf_->getSampleSet();
    const Arrayd&    angles = ss->getAngles0();

    angles0_.insert(angles[0]);

    if (angles.size() == 1) return;

    angles0_.insert(angles[angles.size() - 1]);

    for (int index = 1, prevIndex = 0; index < angles.size() - 1; ++index) {
        bool inserted = false;

        for (int i1 = 0; i1 < ss->getNumAngles1() && !inserted; ++i1) {
        for (int i2 = 0; i2 < ss->getNumAngles2() && !inserted; ++i2) {
        for (int i3 = 0; i3 < ss->getNumAngles3() && !inserted; ++i3) {
            if (hasDownwardDir(*brdf_, index, i1, i2, i3)) break;

            int nextIndex = index + 1;

            double angle = angles[index];
            double prevAngle = angles[prevIndex];
            double nextAngle = angles[nextIndex];

            bool prevDownward = hasDownwardDir(*brdf_, prevIndex, i1, i2, i3);
            bool nextDownward = hasDownwardDir(*brdf_, nextIndex, i1, i2, i3);

            Arrayd sp     = ss->getSpectrum(index,     i1, i2, i3).cast<double>();
            Arrayd prevSp = ss->getSpectrum(prevIndex, i1, i2, i3).cast<double>();
            Arrayd nextSp = ss->getSpectrum(nextIndex, i1, i2, i3).cast<double>();
            
            if (isExtraAngle(angle, prevAngle, nextAngle,
                             prevDownward, nextDownward,
                             sp, prevSp, nextSp)) {
                continue;
            }

            angles0_.insert(angle);
            prevIndex = index;
            inserted = true;
        }}}
    }
}

void Optimizer::setupAngles1()
{
    const SampleSet* ss = brdf_->getSampleSet();
    const Arrayd&    angles = ss->getAngles1();

    angles1_.insert(angles[0]);

    if (angles.size() == 1) return;

    angles1_.insert(angles[angles.size() - 1]);

    for (int index = 1, prevIndex = 0; index < angles.size() - 1; ++index) {
        bool inserted = false;

        for (int i0 = 0; i0 < ss->getNumAngles0() && !inserted; ++i0) {
        for (int i2 = 0; i2 < ss->getNumAngles2() && !inserted; ++i2) {
        for (int i3 = 0; i3 < ss->getNumAngles3() && !inserted; ++i3) {
            if (hasDownwardDir(*brdf_, i0, index, i2, i3)) break;

            int nextIndex = index + 1;

            double angle = angles[index];
            double prevAngle = angles[prevIndex];
            double nextAngle = angles[nextIndex];

            bool prevDownward = hasDownwardDir(*brdf_, i0, prevIndex, i2, i3);
            bool nextDownward = hasDownwardDir(*brdf_, i0, nextIndex, i2, i3);

            Arrayd sp     = ss->getSpectrum(i0, index,     i2, i3).cast<double>();
            Arrayd prevSp = ss->getSpectrum(i0, prevIndex, i2, i3).cast<double>();
            Arrayd nextSp = ss->getSpectrum(i0, nextIndex, i2, i3).cast<double>();
            
            if (isExtraAngle(angle, prevAngle, nextAngle,
                             prevDownward, nextDownward,
                             sp, prevSp, nextSp)) {
                continue;
            }

            angles1_.insert(angle);
            prevIndex = index;
            inserted = true;
        }}}
    }
}

void Optimizer::setupAngles2()
{
    const SampleSet* ss = brdf_->getSampleSet();
    const Arrayd&    angles = ss->getAngles2();

    angles2_.insert(angles[0]);

    if (angles.size() == 1) return;

    angles2_.insert(angles[angles.size() - 1]);

    for (int index = 1, prevIndex = 0; index < angles.size() - 1; ++index) {
        bool inserted = false;

        for (int i0 = 0; i0 < ss->getNumAngles0() && !inserted; ++i0) {
        for (int i1 = 0; i1 < ss->getNumAngles1() && !inserted; ++i1) {
        for (int i3 = 0; i3 < ss->getNumAngles3() && !inserted; ++i3) {
            if (hasDownwardDir(*brdf_, i0, i1, index, i3)) break;

            int nextIndex = index + 1;

            double angle = angles[index];
            double prevAngle = angles[prevIndex];
            double nextAngle = angles[nextIndex];

            bool prevDownward = hasDownwardDir(*brdf_, i0, i1, prevIndex, i3);
            bool nextDownward = hasDownwardDir(*brdf_, i0, i1, nextIndex, i3);

            Arrayd sp     = ss->getSpectrum(i0, i1, index,     i3).cast<double>();
            Arrayd prevSp = ss->getSpectrum(i0, i1, prevIndex, i3).cast<double>();
            Arrayd nextSp = ss->getSpectrum(i0, i1, nextIndex, i3).cast<double>();
            
            if (isExtraAngle(angle, prevAngle, nextAngle,
                             prevDownward, nextDownward,
                             sp, prevSp, nextSp)) {
                continue;
            }

            angles2_.insert(angle);
            prevIndex = index;
            inserted = true;
        }}}
    }
}

void Optimizer::setupAngles3()
{
    const SampleSet* ss = brdf_->getSampleSet();
    const Arrayd&    angles = ss->getAngles3();

    angles3_.insert(angles[0]);

    if (angles.size() == 1) return;

    angles3_.insert(angles[angles.size() - 1]);

    for (int index = 1, prevIndex = 0; index < angles.size() - 1; ++index) {
        bool inserted = false;

        for (int i0 = 0; i0 < ss->getNumAngles0() && !inserted; ++i0) {
        for (int i1 = 0; i1 < ss->getNumAngles1() && !inserted; ++i1) {
        for (int i2 = 0; i2 < ss->getNumAngles2() && !inserted; ++i2) {
            if (hasDownwardDir(*brdf_, i0, i1, i2, index)) break;

            int nextIndex = index + 1;

            double angle = angles[index];
            double prevAngle = angles[prevIndex];
            double nextAngle = angles[nextIndex];

            bool prevDownward = hasDownwardDir(*brdf_, i0, i1, i2, prevIndex);
            bool nextDownward = hasDownwardDir(*brdf_, i0, i1, i2, nextIndex);

            Arrayd sp     = ss->getSpectrum(i0, i1, i2, index).cast<double>();
            Arrayd prevSp = ss->getSpectrum(i0, i1, i2, prevIndex).cast<double>();
            Arrayd nextSp = ss->getSpectrum(i0, i1, i2, nextIndex).cast<double>();
            
            if (isExtraAngle(angle, prevAngle, nextAngle,
                             prevDownward, nextDownward,
                             sp, prevSp, nextSp)) {
                continue;
            }

            angles3_.insert(angle);
            prevIndex = index;
            inserted = true;
        }}}
    }
}

bool Optimizer::isExtraAngle(double        angle,
                             double        prevAngle,
                             double        nextAngle,
                             bool          prevDownward,
                             bool          nextDownward,
                             const Arrayd& sp,
                             const Arrayd& prevSp,
                             const Arrayd& nextSp)
{
    Arrayd diffSp;

    if (!prevDownward && !nextDownward) {
        Arrayd interpolatedSp = lerp(prevSp, nextSp, (angle - prevAngle) / (nextAngle - prevAngle));
        diffSp = (sp - interpolatedSp).abs();
    }
    else if (prevDownward) {
        diffSp = (sp - nextSp).abs();
    }
    else if (nextDownward) {
        diffSp = (sp - prevSp).abs();
    }
    else {
        return true;
    }

    Arrayd ratioSp = diffSp.cwiseQuotient(sp.cwiseMax(EPSILON_F));

    if (diffSp.maxCoeff()  <= diffThreshold_ ||
        ratioSp.maxCoeff() <= ratioThreshold_) {
        return true;
    }

    return false;
}

void Optimizer::updateBrdf()
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

    initializeSpectra<LinearInterpolator>(*origBrdf, brdf_);

    delete origBrdf;
}
