// =================================================================== //
// Copyright (C) 2014-2020 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#include <libbsdf/Common/SpectrumUtility.h>

using namespace lb;

Vec3 SpectrumUtility::spectrumToXyz(const Spectrum& spectrum,
                                    const Arrayf&   wavelengths)
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

    const Vec3d neutralXyz(0.95047, 1, 1.08883);
    const Vec3d normalizingConst = NORMALIZING_CONSTANT_CIE_XYZ.cast<Vec3d::Scalar>();
    sumXyz = sumXyz.cwiseQuotient(normalizingConst).cwiseProduct(neutralXyz);

    return Vec3(sumXyz.cast<Vec3::Scalar>());
}

Vec3 SpectrumUtility::spectrumToXyz(const Spectrum& spectrum,
                                    ColorModel      colorModel,
                                    const Arrayf&   wavelengths)
{
    switch (colorModel) {
        case MONOCHROMATIC_MODEL:
            return Vec3(0.95047f, 1.0f, 1.08883f) * spectrum[0];
        case RGB_MODEL:
            return srgbToXyz(toVec3(spectrum));
        case XYZ_MODEL:
            return toVec3(spectrum);
        case SPECTRAL_MODEL:
            return spectrumToXyz(spectrum, wavelengths);
        default:
            lbError << "[SpectrumUtility::spectrumToY] Invalid color model: " << colorModel;
            return Vec3(0.0, 0.0, 0.0);
    }
}

Vec3::Scalar SpectrumUtility::spectrumToY(const Spectrum& spectrum, const Arrayf& wavelengths)
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

    return static_cast<Vec3::Scalar>(sumY / NORMALIZING_CONSTANT_CIE_XYZ[1]);
}

float SpectrumUtility::spectrumToY(const Spectrum&  spectrum,
                                   ColorModel       colorModel,
                                   const Arrayf&    wavelengths)
{
    switch (colorModel) {
        case MONOCHROMATIC_MODEL:
            return spectrum[0];
        case RGB_MODEL:
            return srgbToXyz<Vec3f>(spectrum)[1];
        case XYZ_MODEL:
            return spectrum[1];
        case SPECTRAL_MODEL:
            return static_cast<float>(spectrumToY(spectrum, wavelengths));
        default:
            lbError << "[SpectrumUtility::spectrumToY] Invalid color model: " << colorModel;
            return 0.0f;
    }
}

Vec3 SpectrumUtility::computeXyzNormalizingConstant()
{
    Vec3 prevXyz(CieData::XYZ[0],
                 CieData::XYZ[1],
                 CieData::XYZ[2]);
    prevXyz *= CieData::D65[0];

    Vec3d sumXyz = Vec3d::Zero();
    const float interval = (CieData::maxWavelength - CieData::minWavelength)
                         / (CieData::numWavelengths - 1);

    // Trapezoidal rule
    for (int i = 1; i < CieData::numWavelengths; ++i) {
        Vec3 xyz(CieData::XYZ[i * 3],
                 CieData::XYZ[i * 3 + 1],
                 CieData::XYZ[i * 3 + 2]);
        xyz *= CieData::D65[i];

        Vec3 area = interval * (prevXyz + xyz);
        sumXyz[0] += area[0];
        sumXyz[1] += area[1];
        sumXyz[2] += area[2];

        prevXyz = xyz;
    }
    sumXyz /= 2.0;

    return sumXyz.cast<Vec3::Scalar>();
}

const Vec3 SpectrumUtility::NORMALIZING_CONSTANT_CIE_XYZ(10043.8f, 10567.3f, 11507.4f);
