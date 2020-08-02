// =================================================================== //
// Copyright (C) 2020 Kimura Ryo                                       //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#include <libbsdf/Common/Utility.h>

#include <chrono>
#include <iomanip>

#include <libbsdf/Common/MunsellData.h>

using namespace lb;

Vec3::Scalar lb::computeCiede2000(const Vec3& lab0, const Vec3& lab1)
{
    using std::atan2;
    using std::cos;
    using std::exp;
    using std::fabs;
    using std::pow;
    using std::sin;
    using std::sqrt;

    double C0 = sqrt(lab0[1] * lab0[1] + lab0[2] * lab0[2]);
    double C1 = sqrt(lab1[1] * lab1[1] + lab1[2] * lab1[2]);
    double avgC = (C0 + C1) * 0.5;
    double G = 0.5 * (1.0 - sqrt(pow(avgC, 7.0) / (pow(avgC, 7.0) + pow(25, 7))));
    double a0Prime = (1.0 + G) * lab0[1];
    double a1Prime = (1.0 + G) * lab1[1];
    double C0Prime = sqrt((a0Prime * a0Prime) + (lab0[2] * lab0[2]));
    double C1Prime = sqrt((a1Prime * a1Prime) + (lab1[2] * lab1[2]));
    double h0Prime;
    if (lab0[2] == 0.0 && a0Prime == 0.0) {
        h0Prime = 0.0;
    }
    else {
        h0Prime = atan2(lab0[2], a0Prime);
        if (h0Prime < 0.0) {
            h0Prime += TAU_D;
        }
    }
    double h1Prime;
    if (lab1[2] == 0.0 && a1Prime == 0.0) {
        h1Prime = 0.0;
    }
    else {
        h1Prime = atan2(lab1[2], a1Prime);
        if (h1Prime < 0.0) {
            h1Prime += TAU_D;
        }
    }

    double deltaLPrime = lab1[0] - lab0[0];
    double deltaCPrime = C1Prime - C0Prime;
    double deltahPrime;
    if ((C0Prime * C1Prime) == 0.0) {
        deltahPrime = 0.0;
    }
    else {
        deltahPrime = h1Prime - h0Prime;
        if (deltahPrime < -PI_D) {
            deltahPrime += TAU_D;
        }
        else if (deltahPrime > PI_D) {
            deltahPrime -= TAU_D;
        }
    }
    double deltaHPrime = 2.0 * sqrt(C0Prime * C1Prime) * sin(deltahPrime * 0.5);

    double avgLPrime = (lab0[0] + lab1[0]) * 0.5;
    double avgCPrime = (C0Prime + C1Prime) * 0.5;
    double avghPrime;
    double sumhPrime = h0Prime + h1Prime;
    if ((C0Prime * C1Prime) == 0.0) {
        avghPrime = sumhPrime;
    }
    else {
        if (fabs(h0Prime - h1Prime) <= PI_D) {
            avghPrime = sumhPrime * 0.5;
        }
        else {
            if (sumhPrime < TAU_D) {
                avghPrime = (sumhPrime + TAU_D) * 0.5;
            }
            else {
                avghPrime = (sumhPrime - TAU_D) * 0.5;
            }
        }
    }

    double T = 1.0 - (0.17 * cos(avghPrime - toRadian(30.0))) +
                     (0.24 * cos(2.0 * avghPrime)) +
                     (0.32 * cos((3.0 * avghPrime) + toRadian(6.0))) -
                     (0.20 * cos((4.0 * avghPrime) - toRadian(63.0)));
    double deltaTheta = toRadian(30.0) * exp(-pow((avghPrime - toRadian(275.0)) / toRadian(25.0), 2.0));
    double RC = 2.0 * sqrt(pow(avgCPrime, 7.0) / (pow(avgCPrime, 7.0) + pow(25, 7)));
    double SL = 1.0 + ((0.015 * pow(avgLPrime - 50.0, 2.0)) / sqrt(20 + pow(avgLPrime - 50.0, 2.0)));
    double SC = 1.0 + (0.045 * avgCPrime);
    double SH = 1.0 + (0.015 * avgCPrime * T);
    double RT = (-sin(2.0 * deltaTheta)) * RC;

    constexpr double kL = 1.0, kC = 1.0, kH = 1.0;

    double deltaE = sqrt(pow(deltaLPrime / (kL * SL), 2.0) +
                         pow(deltaCPrime / (kC * SC), 2.0) +
                         pow(deltaHPrime / (kH * SH), 2.0) +
                         (RT * (deltaCPrime / (kC * SC)) * (deltaHPrime / (kH * SH))));

    return deltaE;
}

Vec3 lb::findMunsellProperties(const Vec3& xyz, std::string* hue, float* value, int* chroma)
{
    Vec3::Scalar minDistance = std::numeric_limits<Vec3::Scalar>::max();
    Vec3 lab = xyzToLab(xyz);

    Vec3 nearestXyz;

    float prevV = -1.0f;
    for (int i = 0; i < MunsellData::numColors; ++i) {
        Vec3::Scalar munsellY = MunsellData::xyY[i * 3 + 2] / MunsellData::NORMALIZING_CONSTANT_Y;
        Vec3 xyy(MunsellData::xyY[i * 3],
                 MunsellData::xyY[i * 3 + 1],
                 munsellY);
        Vec3 munsellXyz = xyyToXyz(xyy);
        Vec3 muncellLab = xyzToLab(munsellXyz);

        Vec3::Scalar distance = computeCiede2000(lab, muncellLab);
        if (distance < minDistance) {
            minDistance = distance;
            *hue    = MunsellData::H[i];
            *value  = MunsellData::V[i];
            *chroma = MunsellData::C[i];
            nearestXyz = munsellXyz;
        }

        // Find a neutral color.
        if (prevV != MunsellData::V[i]) {
            Vec3 neutralXyz = Vec3(0.95047, 1, 1.08883) * munsellY;
            Vec3 neutralLab = xyzToLab(neutralXyz);

            distance = computeCiede2000(lab, neutralLab);
            if (distance < minDistance) {
                minDistance = distance;
                *hue = "N";
                *value = MunsellData::V[i];
                *chroma = 0;
                nearestXyz = neutralXyz;
            }

            prevV = MunsellData::V[i];
        }
    }

    // Find the black color.
    if (*chroma == 0) {
        Vec3 blackLab = Vec3::Zero();
        Vec3::Scalar distance = computeCiede2000(lab, blackLab);
        if (distance < minDistance) {
            minDistance = distance;
            *hue = "N";
            *value = 0.0f;
            nearestXyz.setZero();
        }
    }

    return nearestXyz;
}

std::string lb::getDate()
{
    std::stringstream ss;
    std::time_t t = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
#if _MSC_VER
    std::tm tm;
    if (localtime_s(&tm, &t) == 0) {
        ss << std::put_time(&tm, "%F");
    }
#else
    ss << std::put_time(std::localtime(&t), "%F");
#endif

    return ss.str();
}
