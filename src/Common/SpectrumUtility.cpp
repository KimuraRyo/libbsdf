// =================================================================== //
// Copyright (C) 2014-2017 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#include <libbsdf/Common/SpectrumUtility.h>

using namespace lb;

Vec3 SpectrumUtility::computeNormalizingConstant_sRGB()
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

    return xyzToSrgb(sumXyz.cast<Vec3f::Scalar>());
}

const Vec3 SpectrumUtility::NORMALIZING_CONSTANT_SRGB(10566.4f, 10567.4f, 10568.8f);
const Vec3 SpectrumUtility::NORMALIZING_CONSTANT_CIE_XYZ(10043.8f, 10567.3f, 11507.4f);
