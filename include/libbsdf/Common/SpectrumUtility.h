// =================================================================== //
// Copyright (C) 2014-2019 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#ifndef LIBBSDF_SPECTRUM_UTILITY_H
#define LIBBSDF_SPECTRUM_UTILITY_H

#include <libbsdf/Common/CieData.h>
#include <libbsdf/Common/Utility.h>
#include <libbsdf/Common/Vector.h>

namespace lb {

/*!
 * \class   SpectrumUtility
 * \brief   The SpectrumUtility class provides utility functions for a spectrum.
 */
class SpectrumUtility
{
public:
    /*! Converts from a spectrum to sRGB. Negative sRGB values are clamped. */
    template <typename T>
    static Vec3 spectrumToSrgb(const T& spectrum, const T& wavelengths);

    /*! Converts from a spectrum to Y of CIE-XYZ. */
    template <typename T>
    static Vec3::Scalar spectrumToY(const T& spectrum, const T& wavelengths);

    /*! Converts from a spectrum to CIE-XYZ. */
    template <typename T>
    static Vec3 spectrumToXyz(const T& spectrum, const T& wavelengths);

    /*! Converts from a wavelength to sRGB. Negative or more than 1.0 sRGB values are clamped. */
    static Vec3 wavelengthToSrgb(float wavelength);

private:
    /*! Finds the nearest index in the array of wavelengths. */
    static int findNearestIndex(float wavelength);

    /*! Computes constants to normalize sRGB. */
    static Vec3 computeNormalizingConstant_sRGB();

    /*! Precomputed constants to normalize sRGB values. */
    static const Vec3 NORMALIZING_CONSTANT_SRGB;

    /*! Precomputed constants to normalize CIE_XYZ values. */
    static const Vec3 NORMALIZING_CONSTANT_CIE_XYZ;
};

template <typename T>
Vec3 SpectrumUtility::spectrumToSrgb(const T& spectrum, const T& wavelengths)
{
    Vec3 xyz = spectrumToXyz(spectrum, wavelengths);
    Vec3 rgb = xyzToSrgb(xyz);
    rgb = rgb.cwiseMax(0.0);
    return rgb.cwiseQuotient(NORMALIZING_CONSTANT_SRGB);
}

template <typename T>
Vec3::Scalar SpectrumUtility::spectrumToY(const T& spectrum, const T& wavelengths)
{
    assert(spectrum.size() == wavelengths.size());

    float prevWl = wavelengths[0];
    int index0 = findNearestIndex(prevWl);
    Vec3::Scalar prevY = CieData::XYZ[index0 * 3 + 1];
    prevY *= CieData::D65[index0] * spectrum[0];

    Vec3d::Scalar sumY = 0.0;

    // Trapezoidal rule
    for (int i = 1; i < wavelengths.size(); ++i) {
        float wl = wavelengths[i];
        int index = findNearestIndex(wl);
        Vec3::Scalar y = CieData::XYZ[index * 3 + 1];
        y *= CieData::D65[index] * spectrum[i];

        Vec3::Scalar area = (wl - prevWl) * (prevY + y);
        sumY += area;

        prevWl = wl;
        prevY = y;
    }
    sumY /= 2.0;

    return sumY / NORMALIZING_CONSTANT_CIE_XYZ[1];
}

template <typename T>
Vec3 SpectrumUtility::spectrumToXyz(const T& spectrum, const T& wavelengths)
{
    assert(spectrum.size() == wavelengths.size());

    float prevWl = wavelengths[0];
    int index0 = findNearestIndex(prevWl);
    Vec3 prevXyz(CieData::XYZ[index0 * 3],
                 CieData::XYZ[index0 * 3 + 1],
                 CieData::XYZ[index0 * 3 + 2]);
    prevXyz *= CieData::D65[index0] * spectrum[0];

    Vec3d sumXyz = Vec3d::Zero();

    // Trapezoidal rule
    for (int i = 1; i < wavelengths.size(); ++i) {
        float wl = wavelengths[i];
        int index = findNearestIndex(wl);
        Vec3 xyz(CieData::XYZ[index * 3],
                 CieData::XYZ[index * 3 + 1],
                 CieData::XYZ[index * 3 + 2]);
        xyz *= CieData::D65[index] * spectrum[i];

        Vec3 area = (wl - prevWl) * (prevXyz + xyz);

        sumXyz[0] += area[0];
        sumXyz[1] += area[1];
        sumXyz[2] += area[2];

        prevWl = wl;
        prevXyz = xyz;
    }
    sumXyz /= 2.0;

    return Vec3(sumXyz.cast<Vec3::Scalar>());
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
