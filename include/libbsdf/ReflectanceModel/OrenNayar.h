// =================================================================== //
// Copyright (C) 2016 Kimura Ryo                                       //
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

    OrenNayar(const Vec3&   albedo,
              float         roughness)
              : albedo_(albedo),
                roughness_(roughness)
    {
        parameters_.push_back(Parameter("Albedo",       &albedo_));
        parameters_.push_back(Parameter("Roughness",    &roughness_));
    }

    static Vec3 compute(const Vec3& L,
                        const Vec3& V,
                        const Vec3& N,
                        const Vec3& albedo,
                        float       roughness);

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
    Vec3    albedo_;
    float   roughness_;
};

/*
 * Implementation
 */

inline Vec3 OrenNayar::compute(const Vec3& L,
                               const Vec3& V,
                               const Vec3& N,
                               const Vec3& albedo,
                               float       roughness)
{
    using std::abs;
    using std::acos;
    using std::max;
    using std::min;
    using std::pow;
    using std::sin;
    using std::tan;

    float dotLN = L.dot(N);
    float dotVN = V.dot(N);

    if (dotLN <= 0.0f || dotVN <= 0.0f) return Vec3::Zero();

    float cosPhiDiff;
    if (dotLN == 1.0f || dotVN == 1.0f) {
        cosPhiDiff = 0.0f;
    }
    else {
        Vec3 projectedL = (L - N * dotLN).normalized();
        Vec3 projectedV = (V - N * dotVN).normalized();
        cosPhiDiff = projectedL.dot(projectedV);
    }

    float thetaL = acos(dotLN);
    float thetaV = acos(dotVN);
    float alpha = max(thetaL, thetaV);
    float beta = min(thetaL, thetaV);

    float sqR = roughness * roughness;

    float C1 = 1.0f - 0.5f * sqR / (sqR + 0.33f);

    float C2 = 0.45f * sqR / (sqR + 0.09f);
    if (cosPhiDiff >= 0.0f) {
        C2 *= sin(alpha);
    }
    else {
        C2 *= (sin(alpha) - pow(2.0f * beta / PI_F, 3.0f));
    }

    float alphaBetaPi_C3 = (4.0f * alpha * beta) / (PI_F * PI_F);
    float C3 = 0.125f * sqR / (sqR + 0.09f) * alphaBetaPi_C3 * alphaBetaPi_C3;

    Vec3 L1 = albedo / PI_F
            * (C1 +
               cosPhiDiff * C2 * tan(beta) +
               (1.0f - abs(cosPhiDiff)) * C3 * tan((alpha + beta) / 2.0f));

    float betaPi_L2 = 2.0f * beta / PI_F;
    Vec3 L2 = 0.17f * albedo.cwiseProduct(albedo) / PI_F
            * sqR / (sqR + 0.13f)
            * (1.0f - cosPhiDiff * betaPi_L2 * betaPi_L2);
    
    return L1 + L2;
}

} // namespace lb

#endif // LIBBSDF_OREN_NAYAR_H
