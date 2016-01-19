// =================================================================== //
// Copyright (C) 2014-2015 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#include <libbsdf/Brdf/LinearInterpolator.h>

#include <algorithm>

#include <libbsdf/Brdf/SampleSet2D.h>

using namespace lb;

void LinearInterpolator::getSpectrum(const SampleSet&   samples,
                                     float              angle0,
                                     float              angle1,
                                     float              angle2,
                                     float              angle3,
                                     Spectrum*          spectrum)
{
    const Arrayf& angles0 = samples.getAngles0();
    const Arrayf& angles1 = samples.getAngles1();
    const Arrayf& angles2 = samples.getAngles2();
    const Arrayf& angles3 = samples.getAngles3();

    int lIdx0, lIdx1, lIdx2, lIdx3; // index of the lower bound sample point
    int uIdx0, uIdx1, uIdx2, uIdx3; // index of the upper bound sample point
    Vec4 lowerAngles, upperAngles;

    findBounds(angles0, angle0, samples.isEqualIntervalAngles0(), &lIdx0, &uIdx0, &lowerAngles[0], &upperAngles[0]);
    findBounds(angles1, angle1, samples.isEqualIntervalAngles1(), &lIdx1, &uIdx1, &lowerAngles[1], &upperAngles[1]);
    findBounds(angles2, angle2, samples.isEqualIntervalAngles2(), &lIdx2, &uIdx2, &lowerAngles[2], &upperAngles[2]);
    findBounds(angles3, angle3, samples.isEqualIntervalAngles3(), &lIdx3, &uIdx3, &lowerAngles[3], &upperAngles[3]);

    Vec4 angles(angle0, angle1, angle2, angle3);
    Vec4 intervals = (upperAngles - lowerAngles).cwiseMax(EPSILON_F);
    Vec4 weights = (angles - lowerAngles).cwiseQuotient(intervals);

    const Spectrum& sp0000 = samples.getSpectrum(lIdx0, lIdx1, lIdx2, lIdx3);
    const Spectrum& sp0001 = samples.getSpectrum(lIdx0, lIdx1, lIdx2, uIdx3);
    const Spectrum& sp0010 = samples.getSpectrum(lIdx0, lIdx1, uIdx2, lIdx3);
    const Spectrum& sp0011 = samples.getSpectrum(lIdx0, lIdx1, uIdx2, uIdx3);

    const Spectrum& sp0100 = samples.getSpectrum(lIdx0, uIdx1, lIdx2, lIdx3);
    const Spectrum& sp0101 = samples.getSpectrum(lIdx0, uIdx1, lIdx2, uIdx3);
    const Spectrum& sp0110 = samples.getSpectrum(lIdx0, uIdx1, uIdx2, lIdx3);
    const Spectrum& sp0111 = samples.getSpectrum(lIdx0, uIdx1, uIdx2, uIdx3);

    const Spectrum& sp1000 = samples.getSpectrum(uIdx0, lIdx1, lIdx2, lIdx3);
    const Spectrum& sp1001 = samples.getSpectrum(uIdx0, lIdx1, lIdx2, uIdx3);
    const Spectrum& sp1010 = samples.getSpectrum(uIdx0, lIdx1, uIdx2, lIdx3);
    const Spectrum& sp1011 = samples.getSpectrum(uIdx0, lIdx1, uIdx2, uIdx3);

    const Spectrum& sp1100 = samples.getSpectrum(uIdx0, uIdx1, lIdx2, lIdx3);
    const Spectrum& sp1101 = samples.getSpectrum(uIdx0, uIdx1, lIdx2, uIdx3);
    const Spectrum& sp1110 = samples.getSpectrum(uIdx0, uIdx1, uIdx2, lIdx3);
    const Spectrum& sp1111 = samples.getSpectrum(uIdx0, uIdx1, uIdx2, uIdx3);

    Spectrum sp000 = lerp(sp0000, sp0001, weights[3]);
    Spectrum sp001 = lerp(sp0010, sp0011, weights[3]);
    Spectrum sp010 = lerp(sp0100, sp0101, weights[3]);
    Spectrum sp011 = lerp(sp0110, sp0111, weights[3]);
    Spectrum sp100 = lerp(sp1000, sp1001, weights[3]);
    Spectrum sp101 = lerp(sp1010, sp1011, weights[3]);
    Spectrum sp110 = lerp(sp1100, sp1101, weights[3]);
    Spectrum sp111 = lerp(sp1110, sp1111, weights[3]);

    Spectrum sp00 = lerp(sp000, sp001, weights[2]);
    Spectrum sp01 = lerp(sp010, sp011, weights[2]);
    Spectrum sp10 = lerp(sp100, sp101, weights[2]);
    Spectrum sp11 = lerp(sp110, sp111, weights[2]);

    Spectrum sp0 = lerp(sp00, sp01, weights[1]);
    Spectrum sp1 = lerp(sp10, sp11, weights[1]);

    *spectrum = lerp(sp0, sp1, weights[0]);

    assert(spectrum->allFinite());
}

void LinearInterpolator::getSpectrum(const SampleSet&   samples,
                                     float              angle0,
                                     float              angle2,
                                     float              angle3,
                                     Spectrum*          spectrum)
{
    const Arrayf& angles0 = samples.getAngles0();
    const Arrayf& angles2 = samples.getAngles2();
    const Arrayf& angles3 = samples.getAngles3();

    int lIdx0, lIdx2, lIdx3; // index of the lower bound sample point
    int uIdx0, uIdx2, uIdx3; // index of the upper bound sample point
    Vec4 lowerAngles, upperAngles;

    findBounds(angles0, angle0, samples.isEqualIntervalAngles0(), &lIdx0, &uIdx0, &lowerAngles[0], &upperAngles[0]);
    findBounds(angles2, angle2, samples.isEqualIntervalAngles2(), &lIdx2, &uIdx2, &lowerAngles[2], &upperAngles[2]);
    findBounds(angles3, angle3, samples.isEqualIntervalAngles3(), &lIdx3, &uIdx3, &lowerAngles[3], &upperAngles[3]);

    Vec4 angles(angle0, 0.0, angle2, angle3);
    Vec4 intervals = (upperAngles - lowerAngles).cwiseMax(EPSILON_F);
    Vec4 weights = (angles - lowerAngles).cwiseQuotient(intervals);

    const Spectrum& sp0000 = samples.getSpectrum(lIdx0, lIdx2, lIdx3);
    const Spectrum& sp0001 = samples.getSpectrum(lIdx0, lIdx2, uIdx3);
    const Spectrum& sp0010 = samples.getSpectrum(lIdx0, uIdx2, lIdx3);
    const Spectrum& sp0011 = samples.getSpectrum(lIdx0, uIdx2, uIdx3);

    const Spectrum& sp1000 = samples.getSpectrum(uIdx0, lIdx2, lIdx3);
    const Spectrum& sp1001 = samples.getSpectrum(uIdx0, lIdx2, uIdx3);
    const Spectrum& sp1010 = samples.getSpectrum(uIdx0, uIdx2, lIdx3);
    const Spectrum& sp1011 = samples.getSpectrum(uIdx0, uIdx2, uIdx3);

    Spectrum sp000 = lerp(sp0000, sp0001, weights[3]);
    Spectrum sp001 = lerp(sp0010, sp0011, weights[3]);
    Spectrum sp100 = lerp(sp1000, sp1001, weights[3]);
    Spectrum sp101 = lerp(sp1010, sp1011, weights[3]);

    Spectrum sp00 = lerp(sp000, sp001, weights[2]);
    Spectrum sp10 = lerp(sp100, sp101, weights[2]);

    *spectrum = lerp(sp00, sp10, weights[0]);

    assert(spectrum->allFinite());
}

float LinearInterpolator::getValue(const SampleSet& samples,
                                   float            angle0,
                                   float            angle1,
                                   float            angle2,
                                   float            angle3,
                                   int              wavelengthIndex)
{
    const Arrayf& angles0 = samples.getAngles0();
    const Arrayf& angles1 = samples.getAngles1();
    const Arrayf& angles2 = samples.getAngles2();
    const Arrayf& angles3 = samples.getAngles3();

    int lIdx0, lIdx1, lIdx2, lIdx3; // index of the lower bound sample point
    int uIdx0, uIdx1, uIdx2, uIdx3; // index of the upper bound sample point
    Vec4 lowerAngles, upperAngles;

    findBounds(angles0, angle0, samples.isEqualIntervalAngles0(), &lIdx0, &uIdx0, &lowerAngles[0], &upperAngles[0]);
    findBounds(angles1, angle1, samples.isEqualIntervalAngles1(), &lIdx1, &uIdx1, &lowerAngles[1], &upperAngles[1]);
    findBounds(angles2, angle2, samples.isEqualIntervalAngles2(), &lIdx2, &uIdx2, &lowerAngles[2], &upperAngles[2]);
    findBounds(angles3, angle3, samples.isEqualIntervalAngles3(), &lIdx3, &uIdx3, &lowerAngles[3], &upperAngles[3]);

    Vec4 angles(angle0, angle1, angle2, angle3);
    Vec4 intervals = (upperAngles - lowerAngles).cwiseMax(EPSILON_F);
    Vec4 weights = (angles - lowerAngles).cwiseQuotient(intervals);

    float val0000 = samples.getSpectrum(lIdx0, lIdx1, lIdx2, lIdx3)[wavelengthIndex];
    float val0001 = samples.getSpectrum(lIdx0, lIdx1, lIdx2, uIdx3)[wavelengthIndex];
    float val0010 = samples.getSpectrum(lIdx0, lIdx1, uIdx2, lIdx3)[wavelengthIndex];
    float val0011 = samples.getSpectrum(lIdx0, lIdx1, uIdx2, uIdx3)[wavelengthIndex];

    float val0100 = samples.getSpectrum(lIdx0, uIdx1, lIdx2, lIdx3)[wavelengthIndex];
    float val0101 = samples.getSpectrum(lIdx0, uIdx1, lIdx2, uIdx3)[wavelengthIndex];
    float val0110 = samples.getSpectrum(lIdx0, uIdx1, uIdx2, lIdx3)[wavelengthIndex];
    float val0111 = samples.getSpectrum(lIdx0, uIdx1, uIdx2, uIdx3)[wavelengthIndex];

    float val1000 = samples.getSpectrum(uIdx0, lIdx1, lIdx2, lIdx3)[wavelengthIndex];
    float val1001 = samples.getSpectrum(uIdx0, lIdx1, lIdx2, uIdx3)[wavelengthIndex];
    float val1010 = samples.getSpectrum(uIdx0, lIdx1, uIdx2, lIdx3)[wavelengthIndex];
    float val1011 = samples.getSpectrum(uIdx0, lIdx1, uIdx2, uIdx3)[wavelengthIndex];

    float val1100 = samples.getSpectrum(uIdx0, uIdx1, lIdx2, lIdx3)[wavelengthIndex];
    float val1101 = samples.getSpectrum(uIdx0, uIdx1, lIdx2, uIdx3)[wavelengthIndex];
    float val1110 = samples.getSpectrum(uIdx0, uIdx1, uIdx2, lIdx3)[wavelengthIndex];
    float val1111 = samples.getSpectrum(uIdx0, uIdx1, uIdx2, uIdx3)[wavelengthIndex];

    float val000 = lerp(val0000, val0001, weights[3]);
    float val001 = lerp(val0010, val0011, weights[3]);
    float val010 = lerp(val0100, val0101, weights[3]);
    float val011 = lerp(val0110, val0111, weights[3]);
    float val100 = lerp(val1000, val1001, weights[3]);
    float val101 = lerp(val1010, val1011, weights[3]);
    float val110 = lerp(val1100, val1101, weights[3]);
    float val111 = lerp(val1110, val1111, weights[3]);

    float val00 = lerp(val000, val001, weights[2]);
    float val01 = lerp(val010, val011, weights[2]);
    float val10 = lerp(val100, val101, weights[2]);
    float val11 = lerp(val110, val111, weights[2]);

    float val0 = lerp(val00, val01, weights[1]);
    float val1 = lerp(val10, val11, weights[1]);

    float val = lerp(val0, val1, weights[0]);

    assert(!std::isnan(val) && !std::isinf(val));
    return val;
}

float LinearInterpolator::getValue(const SampleSet& samples,
                                   float            angle0,
                                   float            angle2,
                                   float            angle3,
                                   int              wavelengthIndex)
{
    const Arrayf& angles0 = samples.getAngles0();
    const Arrayf& angles2 = samples.getAngles2();
    const Arrayf& angles3 = samples.getAngles3();

    int lIdx0, lIdx2, lIdx3; // index of the lower bound sample point
    int uIdx0, uIdx2, uIdx3; // index of the upper bound sample point
    Vec4 lowerAngles, upperAngles;

    findBounds(angles0, angle0, samples.isEqualIntervalAngles0(), &lIdx0, &uIdx0, &lowerAngles[0], &upperAngles[0]);
    findBounds(angles2, angle2, samples.isEqualIntervalAngles2(), &lIdx2, &uIdx2, &lowerAngles[2], &upperAngles[2]);
    findBounds(angles3, angle3, samples.isEqualIntervalAngles3(), &lIdx3, &uIdx3, &lowerAngles[3], &upperAngles[3]);

    Vec4 angles(angle0, 0.0, angle2, angle3);
    Vec4 intervals = (upperAngles - lowerAngles).cwiseMax(EPSILON_F);
    Vec4 weights = (angles - lowerAngles).cwiseQuotient(intervals);

    float val0000 = samples.getSpectrum(lIdx0, lIdx2, lIdx3)[wavelengthIndex];
    float val0001 = samples.getSpectrum(lIdx0, lIdx2, uIdx3)[wavelengthIndex];
    float val0010 = samples.getSpectrum(lIdx0, uIdx2, lIdx3)[wavelengthIndex];
    float val0011 = samples.getSpectrum(lIdx0, uIdx2, uIdx3)[wavelengthIndex];

    float val1000 = samples.getSpectrum(uIdx0, lIdx2, lIdx3)[wavelengthIndex];
    float val1001 = samples.getSpectrum(uIdx0, lIdx2, uIdx3)[wavelengthIndex];
    float val1010 = samples.getSpectrum(uIdx0, uIdx2, lIdx3)[wavelengthIndex];
    float val1011 = samples.getSpectrum(uIdx0, uIdx2, uIdx3)[wavelengthIndex];

    float val000 = lerp(val0000, val0001, weights[3]);
    float val001 = lerp(val0010, val0011, weights[3]);
    float val100 = lerp(val1000, val1001, weights[3]);
    float val101 = lerp(val1010, val1011, weights[3]);

    float val00 = lerp(val000, val001, weights[2]);
    float val10 = lerp(val100, val101, weights[2]);

    float val = lerp(val00, val10, weights[0]);

    assert(!std::isnan(val) && !std::isinf(val));
    return val;
}

void LinearInterpolator::getSpectrum(const SampleSet2D& ss2,
                                     float              theta,
                                     float              phi,
                                     Spectrum*          spectrum)
{
    const Arrayf& thetaArray = ss2.getThetaArray();
    const Arrayf& phiArray = ss2.getPhiArray();

    int lIdx0, lIdx1; // index of the lower bound sample point
    int uIdx0, uIdx1; // index of the upper bound sample point
    float lowerAngle0, lowerAngle1;
    float upperAngle0, upperAngle1;

    findBounds(thetaArray, theta, ss2.isEqualIntervalTheta(), &lIdx0, &uIdx0, &lowerAngle0, &upperAngle0);
    findBounds(phiArray,   phi,   ss2.isEqualIntervalPhi(),   &lIdx1, &uIdx1, &lowerAngle1, &upperAngle1);

    const Spectrum& sp00 = ss2.getSpectrum(lIdx0, lIdx1);
    const Spectrum& sp01 = ss2.getSpectrum(lIdx0, uIdx1);
    const Spectrum& sp10 = ss2.getSpectrum(uIdx0, lIdx1);
    const Spectrum& sp11 = ss2.getSpectrum(uIdx0, uIdx1);

    float interval0 = std::max(upperAngle0 - lowerAngle0, EPSILON_F);
    float interval1 = std::max(upperAngle1 - lowerAngle1, EPSILON_F);

    float weight0 = (theta - lowerAngle0) / interval0;
    float weight1 = (phi   - lowerAngle1) / interval1;

    Spectrum sp0 = lerp(sp00, sp01, weight1);
    Spectrum sp1 = lerp(sp10, sp11, weight1);

    *spectrum = lerp(sp0, sp1, weight0);

    assert(spectrum->allFinite());
}

void LinearInterpolator::getSpectrum(const SampleSet2D& ss2,
                                     float              theta,
                                     Spectrum*          spectrum)
{
    const Arrayf& thetaArray = ss2.getThetaArray();

    int lIdx0; // index of the lower bound sample point
    int uIdx0; // index of the upper bound sample point
    float lowerAngle0;
    float upperAngle0;

    findBounds(thetaArray, theta, ss2.isEqualIntervalTheta(), &lIdx0, &uIdx0, &lowerAngle0, &upperAngle0);

    const Spectrum& sp0 = ss2.getSpectrum(lIdx0);
    const Spectrum& sp1 = ss2.getSpectrum(uIdx0);

    float interval0 = std::max(upperAngle0 - lowerAngle0, EPSILON_F);
    float weight0 = (theta - lowerAngle0) / interval0;

    *spectrum = lerp(sp0, sp1, weight0);

    assert(spectrum->allFinite());
}

void LinearInterpolator::findBounds(const Arrayf&   angles,
                                    float           angle,
                                    bool            equalIntervalAngles,
                                    int*            lowerIndex,
                                    int*            upperIndex,
                                    Vec4::Scalar*   lowerAngle,
                                    Vec4::Scalar*   upperAngle)
{
    if (angles.size() == 1) {
        *lowerIndex = 0;
        *upperIndex = 0;
        *lowerAngle = angles[0];
        *upperAngle = angles[0];

        return;
    }

    int backIndex = static_cast<int>(angles.size() - 1);
    if (equalIntervalAngles) {
        // Calculate lower and upper indices.
        *lowerIndex = static_cast<int>(backIndex * angle / (angles[backIndex]));
        *lowerIndex = std::min(*lowerIndex, backIndex - 1);
        *upperIndex = *lowerIndex + 1;
    }
    else {
        // Find lower and upper indices.
        const float* anglePtr = std::lower_bound(&angles[0], &angles[0] + angles.size(), angle);
        *upperIndex = clamp(static_cast<int>(anglePtr - &angles[0]), 1, backIndex);
        *lowerIndex = *upperIndex - 1;
    }

    *lowerAngle = angles[*lowerIndex];
    *upperAngle = angles[*upperIndex];
}
