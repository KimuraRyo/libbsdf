// =================================================================== //
// Copyright (C) 2016-2023 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#ifndef LIBBSDF_OREN_NAYAR_H
#define LIBBSDF_OREN_NAYAR_H

#include <libbsdf/Common/Global.h>
#include <libbsdf/ReflectanceModel/ReflectanceModel.h>

namespace lb {

/*! Oren-Nayar reflectance model. */
class OrenNayar : public ReflectanceModel
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    OrenNayar(const Vec3& albedo, double roughness) : albedo_(albedo), roughness_(roughness)
    {
        parameters_.push_back(Parameter("Albedo", &albedo_));
        parameters_.push_back(Parameter("Roughness", &roughness_));
    }

    static Vec3
    compute(const Vec3& L, const Vec3& V, const Vec3& N, const Vec3& albedo, double roughness);

    Vec3 getValue(const Vec3& inDir, const Vec3& outDir) const
    {
        const Vec3 N = Vec3(0.0, 0.0, 1.0);
        return compute(inDir, outDir, N, albedo_, roughness_);
    }

    bool isIsotropic() const { return true; }

    std::string getName() const { return "Oren-Nayar"; }

    std::string getDescription() const
    {
        std::string reference("Michael Oren and Shree K. Nayar, \"Generalization of Lambert's reflectance model,\" Computer Graphics (SIGGRAPH '94 Proceedings), pp. 239-246, July 1994.");
        return reference;
    }

private:
    Vec3   albedo_;
    double roughness_;
};

/*
 * Implementation
 */

inline Vec3 OrenNayar::compute(const Vec3& L,
                               const Vec3& V,
                               const Vec3& N,
                               const Vec3& albedo,
                               double      roughness)
{
    using std::abs;
    using std::acos;
    using std::max;
    using std::min;
    using std::pow;
    using std::sin;
    using std::tan;

    double dotLN = L.dot(N);
    double dotVN = V.dot(N);

    if (dotLN <= 0 || dotVN <= 0) {
        return Vec3::Zero();
    }

    double cosPhiDiff;
    if (dotLN == 1 || dotVN == 1) {
        cosPhiDiff = 0;
    }
    else {
        Vec3 projectedL = (L - N * dotLN).normalized();
        Vec3 projectedV = (V - N * dotVN).normalized();
        cosPhiDiff = projectedL.dot(projectedV);
    }

    double thetaL = acos(dotLN);
    double thetaV = acos(dotVN);
    double alpha = max(thetaL, thetaV);
    double beta = min(thetaL, thetaV);

    double sqR = roughness * roughness;

    double C1 = 1 - 0.5 * sqR / (sqR + 0.33);

    double C2 = 0.45 * sqR / (sqR + 0.09);
    if (cosPhiDiff >= 0) {
        C2 *= sin(alpha);
    }
    else {
        C2 *= (sin(alpha) - pow(2 * beta / PI_D, 3));
    }

    double alphaBetaPi_C3 = (4 * alpha * beta) / (PI_D * PI_D);
    double C3 = 0.125 * sqR / (sqR + 0.09) * alphaBetaPi_C3 * alphaBetaPi_C3;

    Vec3 L1 =
        albedo / PI_D *
        (C1 + cosPhiDiff * C2 * tan(beta) + (1 - abs(cosPhiDiff)) * C3 * tan((alpha + beta) / 2));

    double betaPi_L2 = 2 * beta / PI_D;
    Vec3   L2 = 0.17 * albedo.cwiseProduct(albedo) / PI_D * sqR / (sqR + 0.13) *
              (1 - cosPhiDiff * betaPi_L2 * betaPi_L2);

    return L1 + L2;
}

} // namespace lb

#endif // LIBBSDF_OREN_NAYAR_H
