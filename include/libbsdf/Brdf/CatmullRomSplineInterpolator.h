// =================================================================== //
// Copyright (C) 2014-2015 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#ifndef LIBBSDF_CATMULL_ROM_SPLINE_INTERPOLATOR_H
#define LIBBSDF_CATMULL_ROM_SPLINE_INTERPOLATOR_H

#include <libbsdf/Brdf/SampleSet.h>
#include <libbsdf/Common/Vector.h>
#include <libbsdf/Common/Utility.h>

namespace lb {

class SampleSet2D;

/*!
 * \class   CatmullRomSplineInterpolator
 * \brief   The CatmullRomSplineInterpolator class provides the functions for Catmull-Rom spline interpolation.
 *
 * \a angle1 is not used for isotropic BRDFs.
 */
class CatmullRomSplineInterpolator
{
public:
    /*! Gets the interpolated spectrum of sample points at a set of angles. */
    static void getSpectrum(const SampleSet&    samples,
                            float               angle0,
                            float               angle1,
                            float               angle2,
                            float               angle3,
                            Spectrum*           spectrum);

    /*! Gets the interpolated spectrum of sample points at a set of angles. */
    static void getSpectrum(const SampleSet&    samples,
                            float               angle0,
                            float               angle2,
                            float               angle3,
                            Spectrum*           spectrum);

    /*! Gets the interpolated value of sample points at a set of angles and the index of wavelength. */
    static float getValue(const SampleSet&  samples,
                          float             angle0,
                          float             angle1,
                          float             angle2,
                          float             angle3,
                          int               wavelengthIndex);

    /*! Gets the interpolated value of sample points at a set of angles and the index of wavelength. */
    static float getValue(const SampleSet&  samples,
                          float             angle0,
                          float             angle2,
                          float             angle3,
                          int               wavelengthIndex);

    /*! Gets the interpolated spectrum of sample points at a set of angles. */
    static void getSpectrum(const SampleSet2D&  ss2,
                            float               theta,
                            float               inPhi,
                            Spectrum*           spectrum);

    /*! Gets the interpolated spectrum of sample points at a angle. */
    static void getSpectrum(const SampleSet2D&  ss2,
                            float               theta,
                            Spectrum*           spectrum);

private:
    /*! Finds four near indices and angles. */
    static void findBounds(const Arrayf&    positions,
                           float            posAngle,
                           bool             equalIntervalPositions,
                           bool             repeatBounds,
                           int*             pos0Index,
                           int*             pos1Index,
                           int*             pos2Index,
                           int*             pos3Index,
                           float*           pos0Angle,
                           float*           pos1Angle,
                           float*           pos2Angle,
                           float*           pos3Angle);

    /*! Interpolates spectra of 2D sample points. */
    static Spectrum interpolate2D(const SampleSet&  samples,
                                  int               index0,
                                  int               index1,
                                  int               pos0Index2,
                                  int               pos1Index2,
                                  int               pos2Index2,
                                  int               pos3Index2,
                                  int               pos0Index3,
                                  int               pos1Index3,
                                  int               pos2Index3,
                                  int               pos3Index3,
                                  float             pos0Angle2,
                                  float             pos1Angle2,
                                  float             pos2Angle2,
                                  float             pos3Angle2,
                                  float             pos0Angle3,
                                  float             pos1Angle3,
                                  float             pos2Angle3,
                                  float             pos3Angle3,
                                  float             angle2,
                                  float             angle3);

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
                               float            pos0Angle2,
                               float            pos1Angle2,
                               float            pos2Angle2,
                               float            pos3Angle2,
                               float            pos0Angle3,
                               float            pos1Angle3,
                               float            pos2Angle3,
                               float            pos3Angle3,
                               float            angle2,
                               float            angle3,
                               int              wavelengthIndex);
};

} // namespace lb

#endif // LIBBSDF_CATMULL_ROM_SPLINE_INTERPOLATOR_H
