// =================================================================== //
// Copyright (C) 2014-2020 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#ifndef LIBBSDF_SPECTRUM_UTILITY_H
#define LIBBSDF_SPECTRUM_UTILITY_H

#include <libbsdf/Common/Array.h>
#include <libbsdf/Common/Utility.h>

namespace lb {

/*!
 * \class   SpectrumUtility
 * \brief   The SpectrumUtility class provides utility functions for a spectrum.
 */
class SpectrumUtility
{
public:
    /*! Converts from a spectrum to CIE XYZ. */
    static Vec3 spectrumToXyz(const Spectrum&   spectrum,
                              const Arrayf&     wavelengths);

    /*! Converts from a spectrum to CIE XYZ. */
    static Vec3 spectrumToXyz(const Spectrum&   spectrum,
                              ColorModel        colorModel,
                              const Arrayf&     wavelengths);

    /*! Converts from a spectrum to sRGB. Negative sRGB values are clamped. */
    static Vec3 spectrumToSrgb(const Spectrum&  spectrum,
                               const Arrayf&    wavelengths);

    /*! Converts from a spectrum to Y of CIE XYZ. */
    static Vec3::Scalar spectrumToY(const Spectrum& spectrum,
                                    const Arrayf&   wavelengths);

    /*! Converts from a spectrum to Y of CIE XYZ using lb::ColorModel. */
    static float spectrumToY(const Spectrum&    spectrum,
                             ColorModel         colorModel,
                             const Arrayf&      wavelengths);

    /*! Converts from a wavelength to sRGB. Negative or more than 1.0 sRGB values are clamped. */
    static Vec3 wavelengthToSrgb(float wavelength);

private:
    /*! Finds the nearest index in the array of wavelengths. */
    static int findNearestIndex(float wavelength);

    /*! Computes constants to normalize CIE XYZ. */
    static Vec3 computeXyzNormalizingConstant();

    /*! Precomputed constants to normalize CIE XYZ values. */
    static const Vec3 NORMALIZING_CONSTANT_CIE_XYZ;
};

inline Vec3 SpectrumUtility::spectrumToSrgb(const Spectrum& spectrum,
                                            const Arrayf&   wavelengths)
{
    Vec3 xyz = spectrumToXyz(spectrum, wavelengths);
    Vec3 rgb = xyzToSrgb(xyz);
    return rgb.cwiseMax(0.0);
}

inline Vec3 SpectrumUtility::wavelengthToSrgb(float wavelength)
{
    int index = findNearestIndex(wavelength);
    Vec3 xyz = Vec3f(CieData::XYZ[index * 3],
                     CieData::XYZ[index * 3 + 1],
                     CieData::XYZ[index * 3 + 2]).cast<Vec3::Scalar>();
    Vec3 rgb = xyzToSrgb(xyz);
    rgb = rgb.cwiseMin(1.0);
    rgb = rgb.cwiseMax(0.0);
    rgb /= std::max(rgb.maxCoeff(), Vec3::Scalar(0.001));
    return rgb;
}

inline int SpectrumUtility::findNearestIndex(float wavelength)
{
    float ratio = (wavelength - CieData::minWavelength)
                / (CieData::maxWavelength - CieData::minWavelength);
    int index = static_cast<int>(CieData::numWavelengths * ratio);
    index = clamp(index, 0, CieData::numWavelengths - 1);
    return index;
}

} // namespace lb

#endif // LIBBSDF_SPECTRUM_UTILITY_H
