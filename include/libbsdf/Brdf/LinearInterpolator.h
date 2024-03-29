// =================================================================== //
// Copyright (C) 2014-2023 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#ifndef LIBBSDF_LINEAR_INTERPOLATOR_H
#define LIBBSDF_LINEAR_INTERPOLATOR_H

#include <libbsdf/Brdf/SampleSet.h>

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
    static Spectrum getSpectrum(const SampleSet2D& ss2, double theta, double phi);

    /*! Gets the interpolated spectrum of sample points at a polar angle. */
    static Spectrum getSpectrum(const SampleSet2D& ss2, double theta);
};

} // namespace lb

#endif // LIBBSDF_LINEAR_INTERPOLATOR_H
