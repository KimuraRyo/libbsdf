// =================================================================== //
// Copyright (C) 2016-2023 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#ifndef LIBBSDF_SIMPLIFIED_OREN_NAYAR_H
#define LIBBSDF_SIMPLIFIED_OREN_NAYAR_H

#include <libbsdf/Common/Global.h>
#include <libbsdf/ReflectanceModel/ReflectanceModel.h>

namespace lb {

/*! Oren-Nayar reflectance model (qualitative model). */
class SimplifiedOrenNayar : public ReflectanceModel
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    SimplifiedOrenNayar(const Vec3& albedo, double roughness)
        : albedo_(albedo), roughness_(roughness)
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

    std::string getName() const { return "Oren-Nayar (qualitative model)"; }

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

inline Vec3 SimplifiedOrenNayar::compute(const Vec3& L,
                                         const Vec3& V,
                                         const Vec3& N,
                                         const Vec3& albedo,
                                         double      roughness)
{
    using std::acos;
    using std::max;
    using std::min;
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

    double inTheta = acos(dotLN);
    double outTheta = acos(dotVN);
    double alpha = max(inTheta, outTheta);
    double beta = min(inTheta, outTheta);

    double sqR = roughness * roughness;

    double A = 1 - 0.5 * sqR / (sqR + 0.33);
    double B = 0.45 * sqR / (sqR + 0.09);

    return albedo / PI_D * (A + B * max(0.0, cosPhiDiff) * sin(alpha) * tan(beta));
}

} // namespace lb

#endif // LIBBSDF_SIMPLIFIED_OREN_NAYAR_H
