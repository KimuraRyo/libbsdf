// =================================================================== //
// Copyright (C) 2014-2020 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#ifndef LIBBSDF_SAMPLER_H
#define LIBBSDF_SAMPLER_H

#include <libbsdf/Brdf/Brdf.h>
#include <libbsdf/Brdf/SampleSet2D.h>

namespace lb {

/*!
 * \class   Sampler
 * \brief   The Sampler class provides sampling functions using interpolation and extrapolation.
 */
class Sampler
{
public:
    /*! Gets the interpolated spectrum of sample points at incoming and outgoing directions. */
    template <typename CoordSysT, typename InterpolatorT>
    static Spectrum getSpectrum(const SampleSet&    samples,
                                const Vec3&         inDir,
                                const Vec3&         outDir);

    /*!
     * Gets the interpolated value of sample points at incoming and outgoing directions
     * and the index of wavelength.
     */
    template <typename CoordSysT, typename InterpolatorT>
    static float getValue(const SampleSet&  samples,
                          const Vec3&       inDir,
                          const Vec3&       outDir,
                          int               wavelengthIndex);

    /*! Gets the interpolated spectrum of sample points at incoming and outgoing directions. */
    template <typename InterpolatorT>
    static Spectrum getSpectrum(const Brdf& brdf,
                                const Vec3& inDir,
                                const Vec3& outDir);

    /*!
     * Gets the interpolated value of sample points at incoming and outgoing directions
     * and the index of wavelength.
     */
    template <typename InterpolatorT>
    static float getValue(const Brdf&   brdf,
                          const Vec3&   inDir,
                          const Vec3&   outDir,
                          int           wavelengthIndex);

    /*! Gets the interpolated spectrum of sample points at an incoming direction. */
    template <typename InterpolatorT>
    static Spectrum getSpectrum(const SampleSet2D&  ss2,
                                const Vec3&         inDir);
};

template <typename CoordSysT, typename InterpolatorT>
inline Spectrum Sampler::getSpectrum(const SampleSet&   samples,
                                     const Vec3&        inDir,
                                     const Vec3&        outDir)
{
    float angle0, angle1, angle2, angle3;
    if (samples.isIsotropic()) {
        CoordSysT::fromXyz(inDir, outDir, &angle0, &angle2, &angle3);
        return InterpolatorT::getSpectrum(samples, angle0, angle2, angle3);
    }
    else {
        CoordSysT::fromXyz(inDir, outDir, &angle0, &angle1, &angle2, &angle3);
        return InterpolatorT::getSpectrum(samples, angle0, angle1, angle2, angle3);
    }
}

template <typename CoordSysT, typename InterpolatorT>
inline float Sampler::getValue(const SampleSet& samples,
                               const Vec3&      inDir,
                               const Vec3&      outDir,
                               int              wavelengthIndex)
{
    float angle0, angle1, angle2, angle3;
    if (samples.isIsotropic()) {
        CoordSysT::fromXyz(inDir, outDir, &angle0, &angle2, &angle3);
        return InterpolatorT::getValue(samples, angle0, angle2, angle3, wavelengthIndex);
    }
    else {
        CoordSysT::fromXyz(inDir, outDir, &angle0, &angle1, &angle2, &angle3);
        return InterpolatorT::getValue(samples, angle0, angle1, angle2, angle3, wavelengthIndex);
    }
}

template <typename InterpolatorT>
inline Spectrum Sampler::getSpectrum(const Brdf&    brdf,
                                     const Vec3&    inDir,
                                     const Vec3&    outDir)
{
    const SampleSet* ss = brdf.getSampleSet();

    float angle0, angle1, angle2, angle3;
    if (ss->isIsotropic()) {
        brdf.fromXyz(inDir, outDir, &angle0, &angle2, &angle3);
        return InterpolatorT::getSpectrum(*ss, angle0, angle2, angle3);
    }
    else {
        brdf.fromXyz(inDir, outDir, &angle0, &angle1, &angle2, &angle3);
        return InterpolatorT::getSpectrum(*ss, angle0, angle1, angle2, angle3);
    }
}

template <typename InterpolatorT>
inline float Sampler::getValue(const Brdf&  brdf,
                               const Vec3&  inDir,
                               const Vec3&  outDir,
                               int          wavelengthIndex)
{
    const SampleSet* ss = brdf.getSampleSet();

    float angle0, angle1, angle2, angle3;
    if (ss->isIsotropic()) {
        brdf.fromXyz(inDir, outDir, &angle0, &angle2, &angle3);
        return InterpolatorT::getValue(*ss, angle0, angle2, angle3, wavelengthIndex);
    }
    else {
        brdf.fromXyz(inDir, outDir, &angle0, &angle1, &angle2, &angle3);
        return InterpolatorT::getValue(*ss, angle0, angle1, angle2, angle3, wavelengthIndex);
    }
}

template <typename InterpolatorT>
inline Spectrum Sampler::getSpectrum(const SampleSet2D& ss2,
                                     const Vec3&        inDir)
{
    float inTheta, inPhi;
    if (ss2.isIsotropic()) {
        inTheta = static_cast<float>(std::acos(inDir[2]));
        return InterpolatorT::getSpectrum(ss2, inTheta);
    }
    else {
        SphericalCoordinateSystem::fromXyz(inDir, &inTheta, &inPhi);
        return InterpolatorT::getSpectrum(ss2, inTheta, inPhi);
    }
}

} // namespace lb

#endif // LIBBSDF_SAMPLER_H
