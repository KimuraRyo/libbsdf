// =================================================================== //
// Copyright (C) 2014-2023 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#ifndef LIBBSDF_SMOOTH_INTERPOLATOR_H
#define LIBBSDF_SMOOTH_INTERPOLATOR_H

#include <libbsdf/Brdf/SampleSet.h>
#include <libbsdf/Brdf/SampleSet2D.h>

#include <libbsdf/Common/CentripetalCatmullRomSpline.h>
#include <libbsdf/Common/MonotoneCubicInterpolation.h>

namespace lb {

/*!
 * \class   SmoothInterpolator
 * \brief   The SmoothInterpolator class provides the functions for smooth interpolation
 *          using four sample points for each angle array.
 *
 * \a angle1 is not used for isotropic BRDFs.
 */
template <typename T>
class SmoothInterpolator
{
public:
    /*! Gets the interpolated spectrum of sample points at a set of angles. */
    static Spectrum getSpectrum(const SampleSet& samples,
                                double           angle0,
                                double           angle1,
                                double           angle2,
                                double           angle3);

    /*! Gets the interpolated spectrum of sample points at a set of angles. */
    static Spectrum
    getSpectrum(const SampleSet& samples, double angle0, double angle2, double angle3);

    /*! Gets the interpolated value of sample points at a set of angles and the index of wavelength. */
    static float getValue(const SampleSet& samples,
                          double           angle0,
                          double           angle1,
                          double           angle2,
                          double           angle3,
                          int              wavelengthIndex);

    /*! Gets the interpolated value of sample points at a set of angles and the index of wavelength. */
    static float getValue(const SampleSet& samples,
                          double           angle0,
                          double           angle2,
                          double           angle3,
                          int              wavelengthIndex);

    /*! Gets the interpolated spectrum of sample points at a set of angles. */
    static Spectrum getSpectrum(const SampleSet2D& ss2, double theta, double inPhi);

    /*! Gets the interpolated spectrum of sample points at a angle. */
    static Spectrum getSpectrum(const SampleSet2D& ss2, double theta);

private:
    /*! Finds four near indices and angles. */
    static void findBounds(const Arrayd& positions,
                           double        posAngle,
                           bool          equalIntervalPositions,
                           bool          repeatBounds,
                           int*          pos0Index,
                           int*          pos1Index,
                           int*          pos2Index,
                           int*          pos3Index,
                           double*       pos0Angle,
                           double*       pos1Angle,
                           double*       pos2Angle,
                           double*       pos3Angle);

    /*! Interpolates spectra of 2D sample points. */
    static Spectrum interpolate2D(const SampleSet& samples,
                                  int              index0,
                                  int              index1,
                                  int              pos0Index2,
                                  int              pos1Index2,
                                  int              pos2Index2,
                                  int              pos3Index2,
                                  int              pos0Index3,
                                  int              pos1Index3,
                                  int              pos2Index3,
                                  int              pos3Index3,
                                  double           pos0Angle2,
                                  double           pos1Angle2,
                                  double           pos2Angle2,
                                  double           pos3Angle2,
                                  double           pos0Angle3,
                                  double           pos1Angle3,
                                  double           pos2Angle3,
                                  double           pos3Angle3,
                                  double           angle2,
                                  double           angle3);

    /*! Interpolates values of 2D sample points. */
    static float interpolate2D(const SampleSet& samples,
                               int              index0,
                               int              index1,
                               int              pos0Index2,
                               int              pos1Index2,
                               int              pos2Index2,
                               int              pos3Index2,
                               int              pos0Index3,
                               int              pos1Index3,
                               int              pos2Index3,
                               int              pos3Index3,
                               double           pos0Angle2,
                               double           pos1Angle2,
                               double           pos2Angle2,
                               double           pos3Angle2,
                               double           pos0Angle3,
                               double           pos1Angle3,
                               double           pos2Angle3,
                               double           pos3Angle3,
                               double           angle2,
                               double           angle3,
                               int              wavelengthIndex);
};

using CatmullRomSplineInterpolator  = SmoothInterpolator<CentripetalCatmullRomSpline>;
using MonotoneCubicInterpolator     = SmoothInterpolator<MonotoneCubicInterpolation>;

template <typename T>
Spectrum SmoothInterpolator<T>::getSpectrum(const SampleSet& samples,
                                            double           angle0,
                                            double           angle1,
                                            double           angle2,
                                            double           angle3)
{
    const Arrayd& angles0 = samples.getAngles0();
    const Arrayd& angles1 = samples.getAngles1();
    const Arrayd& angles2 = samples.getAngles2();
    const Arrayd& angles3 = samples.getAngles3();

    int pos0Idx0, pos0Idx1, pos0Idx2, pos0Idx3;
    int pos1Idx0, pos1Idx1, pos1Idx2, pos1Idx3;
    int pos2Idx0, pos2Idx1, pos2Idx2, pos2Idx3;
    int pos3Idx0, pos3Idx1, pos3Idx2, pos3Idx3;

    double pos0Angle0, pos0Angle1, pos0Angle2, pos0Angle3;
    double pos1Angle0, pos1Angle1, pos1Angle2, pos1Angle3;
    double pos2Angle0, pos2Angle1, pos2Angle2, pos2Angle3;
    double pos3Angle0, pos3Angle1, pos3Angle2, pos3Angle3;

    findBounds(angles0, angle0, samples.isEqualIntervalAngles0(), false,
               &pos0Idx0, &pos1Idx0, &pos2Idx0, &pos3Idx0,
               &pos0Angle0, &pos1Angle0, &pos2Angle0, &pos3Angle0);

    findBounds(angles1, angle1, samples.isEqualIntervalAngles1(), true,
               &pos0Idx1, &pos1Idx1, &pos2Idx1, &pos3Idx1,
               &pos0Angle1, &pos1Angle1, &pos2Angle1, &pos3Angle1);

    findBounds(angles2, angle2, samples.isEqualIntervalAngles2(), false,
               &pos0Idx2, &pos1Idx2, &pos2Idx2, &pos3Idx2,
               &pos0Angle2, &pos1Angle2, &pos2Angle2, &pos3Angle2);

    findBounds(angles3, angle3, samples.isEqualIntervalAngles3(), true,
               &pos0Idx3, &pos1Idx3, &pos2Idx3, &pos3Idx3,
               &pos0Angle3, &pos1Angle3, &pos2Angle3, &pos3Angle3);

    Spectrum sp00 = interpolate2D(samples, pos0Idx0, pos0Idx1,
                                  pos0Idx2, pos1Idx2, pos2Idx2, pos3Idx2,
                                  pos0Idx3, pos1Idx3, pos2Idx3, pos3Idx3,
                                  pos0Angle2, pos1Angle2, pos2Angle2, pos3Angle2,
                                  pos0Angle3, pos1Angle3, pos2Angle3, pos3Angle3,
                                  angle2, angle3);

    Spectrum sp01 = interpolate2D(samples, pos0Idx0, pos1Idx1,
                                  pos0Idx2, pos1Idx2, pos2Idx2, pos3Idx2,
                                  pos0Idx3, pos1Idx3, pos2Idx3, pos3Idx3,
                                  pos0Angle2, pos1Angle2, pos2Angle2, pos3Angle2,
                                  pos0Angle3, pos1Angle3, pos2Angle3, pos3Angle3,
                                  angle2, angle3);

    Spectrum sp02 = interpolate2D(samples, pos0Idx0, pos2Idx1,
                                  pos0Idx2, pos1Idx2, pos2Idx2, pos3Idx2,
                                  pos0Idx3, pos1Idx3, pos2Idx3, pos3Idx3,
                                  pos0Angle2, pos1Angle2, pos2Angle2, pos3Angle2,
                                  pos0Angle3, pos1Angle3, pos2Angle3, pos3Angle3,
                                  angle2, angle3);

    Spectrum sp03 = interpolate2D(samples, pos0Idx0, pos3Idx1,
                                  pos0Idx2, pos1Idx2, pos2Idx2, pos3Idx2,
                                  pos0Idx3, pos1Idx3, pos2Idx3, pos3Idx3,
                                  pos0Angle2, pos1Angle2, pos2Angle2, pos3Angle2,
                                  pos0Angle3, pos1Angle3, pos2Angle3, pos3Angle3,
                                  angle2, angle3);

    Spectrum sp10 = interpolate2D(samples, pos1Idx0, pos0Idx1,
                                  pos0Idx2, pos1Idx2, pos2Idx2, pos3Idx2,
                                  pos0Idx3, pos1Idx3, pos2Idx3, pos3Idx3,
                                  pos0Angle2, pos1Angle2, pos2Angle2, pos3Angle2,
                                  pos0Angle3, pos1Angle3, pos2Angle3, pos3Angle3,
                                  angle2, angle3);

    Spectrum sp11 = interpolate2D(samples, pos1Idx0, pos1Idx1,
                                  pos0Idx2, pos1Idx2, pos2Idx2, pos3Idx2,
                                  pos0Idx3, pos1Idx3, pos2Idx3, pos3Idx3,
                                  pos0Angle2, pos1Angle2, pos2Angle2, pos3Angle2,
                                  pos0Angle3, pos1Angle3, pos2Angle3, pos3Angle3,
                                  angle2, angle3);

    Spectrum sp12 = interpolate2D(samples, pos1Idx0, pos2Idx1,
                                  pos0Idx2, pos1Idx2, pos2Idx2, pos3Idx2,
                                  pos0Idx3, pos1Idx3, pos2Idx3, pos3Idx3,
                                  pos0Angle2, pos1Angle2, pos2Angle2, pos3Angle2,
                                  pos0Angle3, pos1Angle3, pos2Angle3, pos3Angle3,
                                  angle2, angle3);

    Spectrum sp13 = interpolate2D(samples, pos1Idx0, pos3Idx1,
                                  pos0Idx2, pos1Idx2, pos2Idx2, pos3Idx2,
                                  pos0Idx3, pos1Idx3, pos2Idx3, pos3Idx3,
                                  pos0Angle2, pos1Angle2, pos2Angle2, pos3Angle2,
                                  pos0Angle3, pos1Angle3, pos2Angle3, pos3Angle3,
                                  angle2, angle3);

    Spectrum sp20 = interpolate2D(samples, pos2Idx0, pos0Idx1,
                                  pos0Idx2, pos1Idx2, pos2Idx2, pos3Idx2,
                                  pos0Idx3, pos1Idx3, pos2Idx3, pos3Idx3,
                                  pos0Angle2, pos1Angle2, pos2Angle2, pos3Angle2,
                                  pos0Angle3, pos1Angle3, pos2Angle3, pos3Angle3,
                                  angle2, angle3);

    Spectrum sp21 = interpolate2D(samples, pos2Idx0, pos1Idx1,
                                  pos0Idx2, pos1Idx2, pos2Idx2, pos3Idx2,
                                  pos0Idx3, pos1Idx3, pos2Idx3, pos3Idx3,
                                  pos0Angle2, pos1Angle2, pos2Angle2, pos3Angle2,
                                  pos0Angle3, pos1Angle3, pos2Angle3, pos3Angle3,
                                  angle2, angle3);

    Spectrum sp22 = interpolate2D(samples, pos2Idx0, pos2Idx1,
                                  pos0Idx2, pos1Idx2, pos2Idx2, pos3Idx2,
                                  pos0Idx3, pos1Idx3, pos2Idx3, pos3Idx3,
                                  pos0Angle2, pos1Angle2, pos2Angle2, pos3Angle2,
                                  pos0Angle3, pos1Angle3, pos2Angle3, pos3Angle3,
                                  angle2, angle3);

    Spectrum sp23 = interpolate2D(samples, pos2Idx0, pos3Idx1,
                                  pos0Idx2, pos1Idx2, pos2Idx2, pos3Idx2,
                                  pos0Idx3, pos1Idx3, pos2Idx3, pos3Idx3,
                                  pos0Angle2, pos1Angle2, pos2Angle2, pos3Angle2,
                                  pos0Angle3, pos1Angle3, pos2Angle3, pos3Angle3,
                                  angle2, angle3);

    Spectrum sp30 = interpolate2D(samples, pos3Idx0, pos0Idx1,
                                  pos0Idx2, pos1Idx2, pos2Idx2, pos3Idx2,
                                  pos0Idx3, pos1Idx3, pos2Idx3, pos3Idx3,
                                  pos0Angle2, pos1Angle2, pos2Angle2, pos3Angle2,
                                  pos0Angle3, pos1Angle3, pos2Angle3, pos3Angle3,
                                  angle2, angle3);

    Spectrum sp31 = interpolate2D(samples, pos3Idx0, pos1Idx1,
                                  pos0Idx2, pos1Idx2, pos2Idx2, pos3Idx2,
                                  pos0Idx3, pos1Idx3, pos2Idx3, pos3Idx3,
                                  pos0Angle2, pos1Angle2, pos2Angle2, pos3Angle2,
                                  pos0Angle3, pos1Angle3, pos2Angle3, pos3Angle3,
                                  angle2, angle3);

    Spectrum sp32 = interpolate2D(samples, pos3Idx0, pos2Idx1,
                                  pos0Idx2, pos1Idx2, pos2Idx2, pos3Idx2,
                                  pos0Idx3, pos1Idx3, pos2Idx3, pos3Idx3,
                                  pos0Angle2, pos1Angle2, pos2Angle2, pos3Angle2,
                                  pos0Angle3, pos1Angle3, pos2Angle3, pos3Angle3,
                                  angle2, angle3);

    Spectrum sp33 = interpolate2D(samples, pos3Idx0, pos3Idx1,
                                  pos0Idx2, pos1Idx2, pos2Idx2, pos3Idx2,
                                  pos0Idx3, pos1Idx3, pos2Idx3, pos3Idx3,
                                  pos0Angle2, pos1Angle2, pos2Angle2, pos3Angle2,
                                  pos0Angle3, pos1Angle3, pos2Angle3, pos3Angle3,
                                  angle2, angle3);

    Spectrum sp0 = T::interpolate(pos0Angle1, pos1Angle1, pos2Angle1, pos3Angle1, sp00, sp01, sp02, sp03, angle1);
    Spectrum sp1 = T::interpolate(pos0Angle1, pos1Angle1, pos2Angle1, pos3Angle1, sp10, sp11, sp12, sp13, angle1);
    Spectrum sp2 = T::interpolate(pos0Angle1, pos1Angle1, pos2Angle1, pos3Angle1, sp20, sp21, sp22, sp23, angle1);
    Spectrum sp3 = T::interpolate(pos0Angle1, pos1Angle1, pos2Angle1, pos3Angle1, sp30, sp31, sp32, sp33, angle1);

    Spectrum sp = T::interpolate(pos0Angle0, pos1Angle0, pos2Angle0, pos3Angle0, sp0, sp1, sp2, sp3, angle0);
    assert(sp.allFinite());

    return sp;
}

template <typename T>
Spectrum SmoothInterpolator<T>::getSpectrum(const SampleSet& samples,
                                            double           angle0,
                                            double           angle2,
                                            double           angle3)
{
    const Arrayd& angles0 = samples.getAngles0();
    const Arrayd& angles2 = samples.getAngles2();
    const Arrayd& angles3 = samples.getAngles3();

    int pos0Idx0, pos0Idx2, pos0Idx3;
    int pos1Idx0, pos1Idx2, pos1Idx3;
    int pos2Idx0, pos2Idx2, pos2Idx3;
    int pos3Idx0, pos3Idx2, pos3Idx3;

    double pos0Angle0, pos0Angle2, pos0Angle3;
    double pos1Angle0, pos1Angle2, pos1Angle3;
    double pos2Angle0, pos2Angle2, pos2Angle3;
    double pos3Angle0, pos3Angle2, pos3Angle3;

    findBounds(angles0, angle0, samples.isEqualIntervalAngles0(), false,
               &pos0Idx0, &pos1Idx0, &pos2Idx0, &pos3Idx0,
               &pos0Angle0, &pos1Angle0, &pos2Angle0, &pos3Angle0);

    findBounds(angles2, angle2, samples.isEqualIntervalAngles2(), false,
               &pos0Idx2, &pos1Idx2, &pos2Idx2, &pos3Idx2,
               &pos0Angle2, &pos1Angle2, &pos2Angle2, &pos3Angle2);

    findBounds(angles3, angle3, samples.isEqualIntervalAngles3(), true,
               &pos0Idx3, &pos1Idx3, &pos2Idx3, &pos3Idx3,
               &pos0Angle3, &pos1Angle3, &pos2Angle3, &pos3Angle3);

    Spectrum sp0 = interpolate2D(samples, pos0Idx0, 0,
                                 pos0Idx2, pos1Idx2, pos2Idx2, pos3Idx2,
                                 pos0Idx3, pos1Idx3, pos2Idx3, pos3Idx3,
                                 pos0Angle2, pos1Angle2, pos2Angle2, pos3Angle2,
                                 pos0Angle3, pos1Angle3, pos2Angle3, pos3Angle3,
                                 angle2, angle3);

    Spectrum sp1 = interpolate2D(samples, pos1Idx0, 0,
                                 pos0Idx2, pos1Idx2, pos2Idx2, pos3Idx2,
                                 pos0Idx3, pos1Idx3, pos2Idx3, pos3Idx3,
                                 pos0Angle2, pos1Angle2, pos2Angle2, pos3Angle2,
                                 pos0Angle3, pos1Angle3, pos2Angle3, pos3Angle3,
                                 angle2, angle3);

    Spectrum sp2 = interpolate2D(samples, pos2Idx0, 0,
                                 pos0Idx2, pos1Idx2, pos2Idx2, pos3Idx2,
                                 pos0Idx3, pos1Idx3, pos2Idx3, pos3Idx3,
                                 pos0Angle2, pos1Angle2, pos2Angle2, pos3Angle2,
                                 pos0Angle3, pos1Angle3, pos2Angle3, pos3Angle3,
                                 angle2, angle3);

    Spectrum sp3 = interpolate2D(samples, pos3Idx0, 0,
                                 pos0Idx2, pos1Idx2, pos2Idx2, pos3Idx2,
                                 pos0Idx3, pos1Idx3, pos2Idx3, pos3Idx3,
                                 pos0Angle2, pos1Angle2, pos2Angle2, pos3Angle2,
                                 pos0Angle3, pos1Angle3, pos2Angle3, pos3Angle3,
                                 angle2, angle3);

    Spectrum sp = T::interpolate(pos0Angle0, pos1Angle0, pos2Angle0, pos3Angle0, sp0, sp1, sp2, sp3, angle0);
    assert(sp.allFinite());

    return sp;
}

template <typename T>
float SmoothInterpolator<T>::getValue(const SampleSet& samples,
                                      double           angle0,
                                      double           angle1,
                                      double           angle2,
                                      double           angle3,
                                      int              wavelengthIndex)
{
    const Arrayd& angles0 = samples.getAngles0();
    const Arrayd& angles1 = samples.getAngles1();
    const Arrayd& angles2 = samples.getAngles2();
    const Arrayd& angles3 = samples.getAngles3();

    int pos0Idx0, pos0Idx1, pos0Idx2, pos0Idx3;
    int pos1Idx0, pos1Idx1, pos1Idx2, pos1Idx3;
    int pos2Idx0, pos2Idx1, pos2Idx2, pos2Idx3;
    int pos3Idx0, pos3Idx1, pos3Idx2, pos3Idx3;

    double pos0Angle0, pos0Angle1, pos0Angle2, pos0Angle3;
    double pos1Angle0, pos1Angle1, pos1Angle2, pos1Angle3;
    double pos2Angle0, pos2Angle1, pos2Angle2, pos2Angle3;
    double pos3Angle0, pos3Angle1, pos3Angle2, pos3Angle3;

    findBounds(angles0, angle0, samples.isEqualIntervalAngles0(), false,
               &pos0Idx0, &pos1Idx0, &pos2Idx0, &pos3Idx0,
               &pos0Angle0, &pos1Angle0, &pos2Angle0, &pos3Angle0);

    findBounds(angles1, angle1, samples.isEqualIntervalAngles1(), true,
               &pos0Idx1, &pos1Idx1, &pos2Idx1, &pos3Idx1,
               &pos0Angle1, &pos1Angle1, &pos2Angle1, &pos3Angle1);

    findBounds(angles2, angle2, samples.isEqualIntervalAngles2(), false,
               &pos0Idx2, &pos1Idx2, &pos2Idx2, &pos3Idx2,
               &pos0Angle2, &pos1Angle2, &pos2Angle2, &pos3Angle2);

    findBounds(angles3, angle3, samples.isEqualIntervalAngles3(), true,
               &pos0Idx3, &pos1Idx3, &pos2Idx3, &pos3Idx3,
               &pos0Angle3, &pos1Angle3, &pos2Angle3, &pos3Angle3);

    float v00 = interpolate2D(samples, pos0Idx0, pos0Idx1,
                              pos0Idx2, pos1Idx2, pos2Idx2, pos3Idx2,
                              pos0Idx3, pos1Idx3, pos2Idx3, pos3Idx3,
                              pos0Angle2, pos1Angle2, pos2Angle2, pos3Angle2,
                              pos0Angle3, pos1Angle3, pos2Angle3, pos3Angle3,
                              angle2, angle3, wavelengthIndex);

    float v01 = interpolate2D(samples, pos0Idx0, pos1Idx1,
                              pos0Idx2, pos1Idx2, pos2Idx2, pos3Idx2,
                              pos0Idx3, pos1Idx3, pos2Idx3, pos3Idx3,
                              pos0Angle2, pos1Angle2, pos2Angle2, pos3Angle2,
                              pos0Angle3, pos1Angle3, pos2Angle3, pos3Angle3,
                              angle2, angle3, wavelengthIndex);

    float v02 = interpolate2D(samples, pos0Idx0, pos2Idx1,
                              pos0Idx2, pos1Idx2, pos2Idx2, pos3Idx2,
                              pos0Idx3, pos1Idx3, pos2Idx3, pos3Idx3,
                              pos0Angle2, pos1Angle2, pos2Angle2, pos3Angle2,
                              pos0Angle3, pos1Angle3, pos2Angle3, pos3Angle3,
                              angle2, angle3, wavelengthIndex);

    float v03 = interpolate2D(samples, pos0Idx0, pos3Idx1,
                              pos0Idx2, pos1Idx2, pos2Idx2, pos3Idx2,
                              pos0Idx3, pos1Idx3, pos2Idx3, pos3Idx3,
                              pos0Angle2, pos1Angle2, pos2Angle2, pos3Angle2,
                              pos0Angle3, pos1Angle3, pos2Angle3, pos3Angle3,
                              angle2, angle3, wavelengthIndex);

    float v10 = interpolate2D(samples, pos1Idx0, pos0Idx1,
                              pos0Idx2, pos1Idx2, pos2Idx2, pos3Idx2,
                              pos0Idx3, pos1Idx3, pos2Idx3, pos3Idx3,
                              pos0Angle2, pos1Angle2, pos2Angle2, pos3Angle2,
                              pos0Angle3, pos1Angle3, pos2Angle3, pos3Angle3,
                              angle2, angle3, wavelengthIndex);

    float v11 = interpolate2D(samples, pos1Idx0, pos1Idx1,
                              pos0Idx2, pos1Idx2, pos2Idx2, pos3Idx2,
                              pos0Idx3, pos1Idx3, pos2Idx3, pos3Idx3,
                              pos0Angle2, pos1Angle2, pos2Angle2, pos3Angle2,
                              pos0Angle3, pos1Angle3, pos2Angle3, pos3Angle3,
                              angle2, angle3, wavelengthIndex);

    float v12 = interpolate2D(samples, pos1Idx0, pos2Idx1,
                              pos0Idx2, pos1Idx2, pos2Idx2, pos3Idx2,
                              pos0Idx3, pos1Idx3, pos2Idx3, pos3Idx3,
                              pos0Angle2, pos1Angle2, pos2Angle2, pos3Angle2,
                              pos0Angle3, pos1Angle3, pos2Angle3, pos3Angle3,
                              angle2, angle3, wavelengthIndex);

    float v13 = interpolate2D(samples, pos1Idx0, pos3Idx1,
                              pos0Idx2, pos1Idx2, pos2Idx2, pos3Idx2,
                              pos0Idx3, pos1Idx3, pos2Idx3, pos3Idx3,
                              pos0Angle2, pos1Angle2, pos2Angle2, pos3Angle2,
                              pos0Angle3, pos1Angle3, pos2Angle3, pos3Angle3,
                              angle2, angle3, wavelengthIndex);

    float v20 = interpolate2D(samples, pos2Idx0, pos0Idx1,
                              pos0Idx2, pos1Idx2, pos2Idx2, pos3Idx2,
                              pos0Idx3, pos1Idx3, pos2Idx3, pos3Idx3,
                              pos0Angle2, pos1Angle2, pos2Angle2, pos3Angle2,
                              pos0Angle3, pos1Angle3, pos2Angle3, pos3Angle3,
                              angle2, angle3, wavelengthIndex);

    float v21 = interpolate2D(samples, pos2Idx0, pos1Idx1,
                              pos0Idx2, pos1Idx2, pos2Idx2, pos3Idx2,
                              pos0Idx3, pos1Idx3, pos2Idx3, pos3Idx3,
                              pos0Angle2, pos1Angle2, pos2Angle2, pos3Angle2,
                              pos0Angle3, pos1Angle3, pos2Angle3, pos3Angle3,
                              angle2, angle3, wavelengthIndex);

    float v22 = interpolate2D(samples, pos2Idx0, pos2Idx1,
                              pos0Idx2, pos1Idx2, pos2Idx2, pos3Idx2,
                              pos0Idx3, pos1Idx3, pos2Idx3, pos3Idx3,
                              pos0Angle2, pos1Angle2, pos2Angle2, pos3Angle2,
                              pos0Angle3, pos1Angle3, pos2Angle3, pos3Angle3,
                              angle2, angle3, wavelengthIndex);

    float v23 = interpolate2D(samples, pos2Idx0, pos3Idx1,
                              pos0Idx2, pos1Idx2, pos2Idx2, pos3Idx2,
                              pos0Idx3, pos1Idx3, pos2Idx3, pos3Idx3,
                              pos0Angle2, pos1Angle2, pos2Angle2, pos3Angle2,
                              pos0Angle3, pos1Angle3, pos2Angle3, pos3Angle3,
                              angle2, angle3, wavelengthIndex);

    float v30 = interpolate2D(samples, pos3Idx0, pos0Idx1,
                              pos0Idx2, pos1Idx2, pos2Idx2, pos3Idx2,
                              pos0Idx3, pos1Idx3, pos2Idx3, pos3Idx3,
                              pos0Angle2, pos1Angle2, pos2Angle2, pos3Angle2,
                              pos0Angle3, pos1Angle3, pos2Angle3, pos3Angle3,
                              angle2, angle3, wavelengthIndex);

    float v31 = interpolate2D(samples, pos3Idx0, pos1Idx1,
                              pos0Idx2, pos1Idx2, pos2Idx2, pos3Idx2,
                              pos0Idx3, pos1Idx3, pos2Idx3, pos3Idx3,
                              pos0Angle2, pos1Angle2, pos2Angle2, pos3Angle2,
                              pos0Angle3, pos1Angle3, pos2Angle3, pos3Angle3,
                              angle2, angle3, wavelengthIndex);

    float v32 = interpolate2D(samples, pos3Idx0, pos2Idx1,
                              pos0Idx2, pos1Idx2, pos2Idx2, pos3Idx2,
                              pos0Idx3, pos1Idx3, pos2Idx3, pos3Idx3,
                              pos0Angle2, pos1Angle2, pos2Angle2, pos3Angle2,
                              pos0Angle3, pos1Angle3, pos2Angle3, pos3Angle3,
                              angle2, angle3, wavelengthIndex);

    float v33 = interpolate2D(samples, pos3Idx0, pos3Idx1,
                              pos0Idx2, pos1Idx2, pos2Idx2, pos3Idx2,
                              pos0Idx3, pos1Idx3, pos2Idx3, pos3Idx3,
                              pos0Angle2, pos1Angle2, pos2Angle2, pos3Angle2,
                              pos0Angle3, pos1Angle3, pos2Angle3, pos3Angle3,
                              angle2, angle3, wavelengthIndex);

    float v0 = T::interpolate(Vec2(pos0Angle1, v00), Vec2(pos1Angle1, v01), Vec2(pos2Angle1, v02), Vec2(pos3Angle1, v03), angle1);
    float v1 = T::interpolate(Vec2(pos0Angle1, v10), Vec2(pos1Angle1, v11), Vec2(pos2Angle1, v12), Vec2(pos3Angle1, v13), angle1);
    float v2 = T::interpolate(Vec2(pos0Angle1, v20), Vec2(pos1Angle1, v21), Vec2(pos2Angle1, v22), Vec2(pos3Angle1, v23), angle1);
    float v3 = T::interpolate(Vec2(pos0Angle1, v30), Vec2(pos1Angle1, v31), Vec2(pos2Angle1, v32), Vec2(pos3Angle1, v33), angle1);

    float val = T::interpolate(Vec2(pos0Angle0, v0), Vec2(pos1Angle0, v1), Vec2(pos2Angle0, v2), Vec2(pos3Angle0, v3), angle0);

    assert(!std::isnan(val) && !std::isinf(val));
    return val;
}

template <typename T>
float SmoothInterpolator<T>::getValue(const SampleSet& samples,
                                      double           angle0,
                                      double           angle2,
                                      double           angle3,
                                      int              wavelengthIndex)
{
    const Arrayd& angles0 = samples.getAngles0();
    const Arrayd& angles2 = samples.getAngles2();
    const Arrayd& angles3 = samples.getAngles3();

    int pos0Idx0, pos0Idx2, pos0Idx3;
    int pos1Idx0, pos1Idx2, pos1Idx3;
    int pos2Idx0, pos2Idx2, pos2Idx3;
    int pos3Idx0, pos3Idx2, pos3Idx3;

    double pos0Angle0, pos0Angle2, pos0Angle3;
    double pos1Angle0, pos1Angle2, pos1Angle3;
    double pos2Angle0, pos2Angle2, pos2Angle3;
    double pos3Angle0, pos3Angle2, pos3Angle3;

    findBounds(angles0, angle0, samples.isEqualIntervalAngles0(), false,
               &pos0Idx0, &pos1Idx0, &pos2Idx0, &pos3Idx0,
               &pos0Angle0, &pos1Angle0, &pos2Angle0, &pos3Angle0);

    findBounds(angles2, angle2, samples.isEqualIntervalAngles2(), false,
               &pos0Idx2, &pos1Idx2, &pos2Idx2, &pos3Idx2,
               &pos0Angle2, &pos1Angle2, &pos2Angle2, &pos3Angle2);

    findBounds(angles3, angle3, samples.isEqualIntervalAngles3(), true,
               &pos0Idx3, &pos1Idx3, &pos2Idx3, &pos3Idx3,
               &pos0Angle3, &pos1Angle3, &pos2Angle3, &pos3Angle3);

    float v0 = interpolate2D(samples, pos0Idx0, 0,
                             pos0Idx2, pos1Idx2, pos2Idx2, pos3Idx2,
                             pos0Idx3, pos1Idx3, pos2Idx3, pos3Idx3,
                             pos0Angle2, pos1Angle2, pos2Angle2, pos3Angle2,
                             pos0Angle3, pos1Angle3, pos2Angle3, pos3Angle3,
                             angle2, angle3, wavelengthIndex);

    float v1 = interpolate2D(samples, pos1Idx0, 0,
                             pos0Idx2, pos1Idx2, pos2Idx2, pos3Idx2,
                             pos0Idx3, pos1Idx3, pos2Idx3, pos3Idx3,
                             pos0Angle2, pos1Angle2, pos2Angle2, pos3Angle2,
                             pos0Angle3, pos1Angle3, pos2Angle3, pos3Angle3,
                             angle2, angle3, wavelengthIndex);

    float v2 = interpolate2D(samples, pos2Idx0, 0,
                             pos0Idx2, pos1Idx2, pos2Idx2, pos3Idx2,
                             pos0Idx3, pos1Idx3, pos2Idx3, pos3Idx3,
                             pos0Angle2, pos1Angle2, pos2Angle2, pos3Angle2,
                             pos0Angle3, pos1Angle3, pos2Angle3, pos3Angle3,
                             angle2, angle3, wavelengthIndex);

    float v3 = interpolate2D(samples, pos3Idx0, 0,
                             pos0Idx2, pos1Idx2, pos2Idx2, pos3Idx2,
                             pos0Idx3, pos1Idx3, pos2Idx3, pos3Idx3,
                             pos0Angle2, pos1Angle2, pos2Angle2, pos3Angle2,
                             pos0Angle3, pos1Angle3, pos2Angle3, pos3Angle3,
                             angle2, angle3, wavelengthIndex);

    float val = T::interpolate(Vec2(pos0Angle0, v0), Vec2(pos1Angle0, v1), Vec2(pos2Angle0, v2), Vec2(pos3Angle0, v3), angle0);

    assert(!std::isnan(val) && !std::isinf(val));
    return val;
}

template <typename T>
Spectrum SmoothInterpolator<T>::getSpectrum(const SampleSet2D& ss2, double theta, double phi)
{
    const Arrayd& thetaArray = ss2.getThetaArray();
    const Arrayd& phiArray = ss2.getPhiArray();

    int pos0Idx0, pos0Idx1;
    int pos1Idx0, pos1Idx1;
    int pos2Idx0, pos2Idx1;
    int pos3Idx0, pos3Idx1;

    double pos0Angle0, pos0Angle1;
    double pos1Angle0, pos1Angle1;
    double pos2Angle0, pos2Angle1;
    double pos3Angle0, pos3Angle1;

    findBounds(thetaArray, theta, ss2.isEqualIntervalTheta(), false,
               &pos0Idx0, &pos1Idx0, &pos2Idx0, &pos3Idx0,
               &pos0Angle0, &pos1Angle0, &pos2Angle0, &pos3Angle0);

    findBounds(phiArray, phi, ss2.isEqualIntervalPhi(), true,
               &pos0Idx1, &pos1Idx1, &pos2Idx1, &pos3Idx1,
               &pos0Angle1, &pos1Angle1, &pos2Angle1, &pos3Angle1);

    const Spectrum& sp00 = ss2.getSpectrum(pos0Idx0, pos0Idx1);
    const Spectrum& sp01 = ss2.getSpectrum(pos0Idx0, pos1Idx1);
    const Spectrum& sp02 = ss2.getSpectrum(pos0Idx0, pos2Idx1);
    const Spectrum& sp03 = ss2.getSpectrum(pos0Idx0, pos3Idx1);

    const Spectrum& sp10 = ss2.getSpectrum(pos1Idx0, pos0Idx1);
    const Spectrum& sp11 = ss2.getSpectrum(pos1Idx0, pos1Idx1);
    const Spectrum& sp12 = ss2.getSpectrum(pos1Idx0, pos2Idx1);
    const Spectrum& sp13 = ss2.getSpectrum(pos1Idx0, pos3Idx1);

    const Spectrum& sp20 = ss2.getSpectrum(pos2Idx0, pos0Idx1);
    const Spectrum& sp21 = ss2.getSpectrum(pos2Idx0, pos1Idx1);
    const Spectrum& sp22 = ss2.getSpectrum(pos2Idx0, pos2Idx1);
    const Spectrum& sp23 = ss2.getSpectrum(pos2Idx0, pos3Idx1);

    const Spectrum& sp30 = ss2.getSpectrum(pos3Idx0, pos0Idx1);
    const Spectrum& sp31 = ss2.getSpectrum(pos3Idx0, pos1Idx1);
    const Spectrum& sp32 = ss2.getSpectrum(pos3Idx0, pos2Idx1);
    const Spectrum& sp33 = ss2.getSpectrum(pos3Idx0, pos3Idx1);

    Spectrum sp0 = T::interpolate(pos0Angle1, pos1Angle1, pos2Angle1, pos3Angle1, sp00, sp01, sp02, sp03, phi);
    Spectrum sp1 = T::interpolate(pos0Angle1, pos1Angle1, pos2Angle1, pos3Angle1, sp10, sp11, sp12, sp13, phi);
    Spectrum sp2 = T::interpolate(pos0Angle1, pos1Angle1, pos2Angle1, pos3Angle1, sp20, sp21, sp22, sp23, phi);
    Spectrum sp3 = T::interpolate(pos0Angle1, pos1Angle1, pos2Angle1, pos3Angle1, sp30, sp31, sp32, sp33, phi);

    Spectrum sp = T::interpolate(pos0Angle0, pos1Angle0, pos2Angle0, pos3Angle0, sp0, sp1, sp2, sp3, theta);
    assert(sp.allFinite());

    return sp;
}

template <typename T>
Spectrum SmoothInterpolator<T>::getSpectrum(const SampleSet2D& ss2, double theta)
{
    const Arrayd& thetaArray = ss2.getThetaArray();

    int pos0Idx0;
    int pos1Idx0;
    int pos2Idx0;
    int pos3Idx0;

    double pos0Angle0;
    double pos1Angle0;
    double pos2Angle0;
    double pos3Angle0;

    findBounds(thetaArray, theta, ss2.isEqualIntervalTheta(), false,
               &pos0Idx0, &pos1Idx0, &pos2Idx0, &pos3Idx0,
               &pos0Angle0, &pos1Angle0, &pos2Angle0, &pos3Angle0);

    const Spectrum& sp0 = ss2.getSpectrum(pos0Idx0);
    const Spectrum& sp1 = ss2.getSpectrum(pos1Idx0);
    const Spectrum& sp2 = ss2.getSpectrum(pos2Idx0);
    const Spectrum& sp3 = ss2.getSpectrum(pos3Idx0);

    Spectrum sp = T::interpolate(pos0Angle0, pos1Angle0, pos2Angle0, pos3Angle0, sp0, sp1, sp2, sp3, theta);
    assert(sp.allFinite());

    return sp;
}

template <typename T>
void SmoothInterpolator<T>::findBounds(const Arrayd& positions,
                                       double        posAngle,
                                       bool          equalIntervalPositions,
                                       bool          repeatBounds,
                                       int*          pos0Index,
                                       int*          pos1Index,
                                       int*          pos2Index,
                                       int*          pos3Index,
                                       double*       pos0Angle,
                                       double*       pos1Angle,
                                       double*       pos2Angle,
                                       double*       pos3Angle)
{
    using std::min;
    using std::max;

    if (positions.size() == 1) {
        *pos0Index = 0;
        *pos1Index = 0;
        *pos2Index = 0;
        *pos3Index = 0;
        *pos0Angle = positions[0];
        *pos1Angle = positions[0];
        *pos2Angle = positions[0];
        *pos3Angle = positions[0];

        return;
    }

    array_util::findBounds(positions, posAngle, equalIntervalPositions, pos1Index, pos2Index, pos1Angle, pos2Angle);

    int backIndex = static_cast<int>(positions.size() - 1);
    if (repeatBounds) {
        if (*pos1Index == 0) {
            *pos0Index = backIndex - 1;
            *pos0Angle = positions[*pos0Index] - positions[backIndex];
        }
        else {
            *pos0Index = *pos1Index - 1;
            *pos0Angle = positions[*pos0Index];
        }

        if (*pos2Index == backIndex) {
            *pos3Index = 1;
            *pos3Angle = positions[*pos3Index] + positions[backIndex];
        }
        else {
            *pos3Index = *pos2Index + 1;
            *pos3Angle = positions[*pos3Index];
        }
    }
    else {
        *pos0Index = max(*pos1Index - 1, 0);
        *pos3Index = min(*pos2Index + 1, backIndex);

        *pos0Angle = positions[*pos0Index];
        *pos3Angle = positions[*pos3Index];
    }
}

template <typename T>
Spectrum SmoothInterpolator<T>::interpolate2D(const SampleSet& samples,
                                              int              index0,
                                              int              index1,
                                              int              pos0Index2,
                                              int              pos1Index2,
                                              int              pos2Index2,
                                              int              pos3Index2,
                                              int              pos0Index3,
                                              int              pos1Index3,
                                              int              pos2Index3,
                                              int              pos3Index3,
                                              double           pos0Angle2,
                                              double           pos1Angle2,
                                              double           pos2Angle2,
                                              double           pos3Angle2,
                                              double           pos0Angle3,
                                              double           pos1Angle3,
                                              double           pos2Angle3,
                                              double           pos3Angle3,
                                              double           angle2,
                                              double           angle3)
{
    const Spectrum& sp00 = samples.getSpectrum(index0, index1, pos0Index2, pos0Index3);
    const Spectrum& sp01 = samples.getSpectrum(index0, index1, pos0Index2, pos1Index3);
    const Spectrum& sp02 = samples.getSpectrum(index0, index1, pos0Index2, pos2Index3);
    const Spectrum& sp03 = samples.getSpectrum(index0, index1, pos0Index2, pos3Index3);

    const Spectrum& sp10 = samples.getSpectrum(index0, index1, pos1Index2, pos0Index3);
    const Spectrum& sp11 = samples.getSpectrum(index0, index1, pos1Index2, pos1Index3);
    const Spectrum& sp12 = samples.getSpectrum(index0, index1, pos1Index2, pos2Index3);
    const Spectrum& sp13 = samples.getSpectrum(index0, index1, pos1Index2, pos3Index3);

    const Spectrum& sp20 = samples.getSpectrum(index0, index1, pos2Index2, pos0Index3);
    const Spectrum& sp21 = samples.getSpectrum(index0, index1, pos2Index2, pos1Index3);
    const Spectrum& sp22 = samples.getSpectrum(index0, index1, pos2Index2, pos2Index3);
    const Spectrum& sp23 = samples.getSpectrum(index0, index1, pos2Index2, pos3Index3);

    const Spectrum& sp30 = samples.getSpectrum(index0, index1, pos3Index2, pos0Index3);
    const Spectrum& sp31 = samples.getSpectrum(index0, index1, pos3Index2, pos1Index3);
    const Spectrum& sp32 = samples.getSpectrum(index0, index1, pos3Index2, pos2Index3);
    const Spectrum& sp33 = samples.getSpectrum(index0, index1, pos3Index2, pos3Index3);

    Spectrum sp0 = T::interpolate(pos0Angle3, pos1Angle3, pos2Angle3, pos3Angle3, sp00, sp01, sp02, sp03, angle3);
    Spectrum sp1 = T::interpolate(pos0Angle3, pos1Angle3, pos2Angle3, pos3Angle3, sp10, sp11, sp12, sp13, angle3);
    Spectrum sp2 = T::interpolate(pos0Angle3, pos1Angle3, pos2Angle3, pos3Angle3, sp20, sp21, sp22, sp23, angle3);
    Spectrum sp3 = T::interpolate(pos0Angle3, pos1Angle3, pos2Angle3, pos3Angle3, sp30, sp31, sp32, sp33, angle3);

    Spectrum sp = T::interpolate(pos0Angle2, pos1Angle2, pos2Angle2, pos3Angle2, sp0, sp1, sp2, sp3, angle2);
    assert(sp.allFinite());

    return sp;
}

template <typename T>
float SmoothInterpolator<T>::interpolate2D(const SampleSet& samples,
                                           int              index0,
                                           int              index1,
                                           int              pos0Index2,
                                           int              pos1Index2,
                                           int              pos2Index2,
                                           int              pos3Index2,
                                           int              pos0Index3,
                                           int              pos1Index3,
                                           int              pos2Index3,
                                           int              pos3Index3,
                                           double           pos0Angle2,
                                           double           pos1Angle2,
                                           double           pos2Angle2,
                                           double           pos3Angle2,
                                           double           pos0Angle3,
                                           double           pos1Angle3,
                                           double           pos2Angle3,
                                           double           pos3Angle3,
                                           double           angle2,
                                           double           angle3,
                                           int              wavelengthIndex)
{
    float v00 = samples.getSpectrum(index0, index1, pos0Index2, pos0Index3)[wavelengthIndex];
    float v01 = samples.getSpectrum(index0, index1, pos0Index2, pos1Index3)[wavelengthIndex];
    float v02 = samples.getSpectrum(index0, index1, pos0Index2, pos2Index3)[wavelengthIndex];
    float v03 = samples.getSpectrum(index0, index1, pos0Index2, pos3Index3)[wavelengthIndex];

    float v10 = samples.getSpectrum(index0, index1, pos1Index2, pos0Index3)[wavelengthIndex];
    float v11 = samples.getSpectrum(index0, index1, pos1Index2, pos1Index3)[wavelengthIndex];
    float v12 = samples.getSpectrum(index0, index1, pos1Index2, pos2Index3)[wavelengthIndex];
    float v13 = samples.getSpectrum(index0, index1, pos1Index2, pos3Index3)[wavelengthIndex];

    float v20 = samples.getSpectrum(index0, index1, pos2Index2, pos0Index3)[wavelengthIndex];
    float v21 = samples.getSpectrum(index0, index1, pos2Index2, pos1Index3)[wavelengthIndex];
    float v22 = samples.getSpectrum(index0, index1, pos2Index2, pos2Index3)[wavelengthIndex];
    float v23 = samples.getSpectrum(index0, index1, pos2Index2, pos3Index3)[wavelengthIndex];

    float v30 = samples.getSpectrum(index0, index1, pos3Index2, pos0Index3)[wavelengthIndex];
    float v31 = samples.getSpectrum(index0, index1, pos3Index2, pos1Index3)[wavelengthIndex];
    float v32 = samples.getSpectrum(index0, index1, pos3Index2, pos2Index3)[wavelengthIndex];
    float v33 = samples.getSpectrum(index0, index1, pos3Index2, pos3Index3)[wavelengthIndex];

    float v0 = T::interpolate(Vec2(pos0Angle3, v00), Vec2(pos1Angle3, v01), Vec2(pos2Angle3, v02), Vec2(pos3Angle3, v03), angle3);
    float v1 = T::interpolate(Vec2(pos0Angle3, v10), Vec2(pos1Angle3, v11), Vec2(pos2Angle3, v12), Vec2(pos3Angle3, v13), angle3);
    float v2 = T::interpolate(Vec2(pos0Angle3, v20), Vec2(pos1Angle3, v21), Vec2(pos2Angle3, v22), Vec2(pos3Angle3, v23), angle3);
    float v3 = T::interpolate(Vec2(pos0Angle3, v30), Vec2(pos1Angle3, v31), Vec2(pos2Angle3, v32), Vec2(pos3Angle3, v33), angle3);

    return T::interpolate(Vec2(pos0Angle2, v0), Vec2(pos1Angle2, v1), Vec2(pos2Angle2, v2), Vec2(pos3Angle2, v3), angle2);
}

} // namespace lb

#endif // LIBBSDF_SMOOTH_INTERPOLATOR_H
