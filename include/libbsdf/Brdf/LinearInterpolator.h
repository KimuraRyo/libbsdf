// =================================================================== //
// Copyright (C) 2014-2015 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#ifndef LIBBSDF_LINEAR_INTERPOLATOR_H
#define LIBBSDF_LINEAR_INTERPOLATOR_H

#include <libbsdf/Brdf/SampleSet.h>
#include <libbsdf/Common/Vector.h>
#include <libbsdf/Common/Utility.h>

namespace lb {

class SampleSet2D;

/*!
 * \class   LinearInterpolator
 * \brief   The LinearInterpolator class provides the functions for linear interpolation.
 *
 * \a angle1 is not used for isotropic BRDFs.
 */
class LinearInterpolator
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
                            float               phi,
                            Spectrum*           spectrum);

    /*! Gets the interpolated spectrum of sample points at a polar angle. */
    static void getSpectrum(const SampleSet2D&  ss2,
                            float               theta,
                            Spectrum*           spectrum);

    /*!
     * Finds neighbor indices and angles.
     *
     * \param lowerIndex Found index of the sample point at the lower bound.
     * \param upperIndex Found index of the sample point at the upper bound.
     * \param lowerAngle Found angle of the sample point at the lower bound.
     * \param upperAngle Found angle of the sample point at the upper bound.
     */
    static void findBounds(const Arrayf&    angles,
                           float            angle,
                           bool             equalIntervalAngles,
                           int*             lowerIndex,
                           int*             upperIndex,
                           Vec4::Scalar*    lowerAngle,
                           Vec4::Scalar*    upperAngle);
};

} // namespace lb

#endif // LIBBSDF_LINEAR_INTERPOLATOR_H
