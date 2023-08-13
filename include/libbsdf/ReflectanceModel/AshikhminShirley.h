// =================================================================== //
// Copyright (C) 2016-2023 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#ifndef LIBBSDF_ASHIKHMIN_SHIRLEY_H
#define LIBBSDF_ASHIKHMIN_SHIRLEY_H

#include <libbsdf/Common/Global.h>
#include <libbsdf/ReflectanceModel/Fresnel.h>
#include <libbsdf/ReflectanceModel/ReflectanceModel.h>

namespace lb {

/*! Ashikhmin-Shirley reflectance model. */
class AshikhminShirley : public ReflectanceModel
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    AshikhminShirley(const Vec3& specularColor,
                     const Vec3& diffuseColor,
                     double      shininessX,
                     double      shininessY)
        : specularColor_(specularColor),
          diffuseColor_(diffuseColor),
          shininessX_(shininessX),
          shininessY_(shininessY)
    {
        parameters_.push_back(Parameter("Specular color", &specularColor_));
        parameters_.push_back(Parameter("Shininess X", &shininessX_, 0.0, 1000.0));
        parameters_.push_back(Parameter("Shininess Y", &shininessY_, 0.0, 1000.0));
        parameters_.push_back(Parameter("Diffuse color", &diffuseColor_));
    }

    static Vec3 compute(const Vec3& L,
                        const Vec3& V,
                        const Vec3& N,
                        const Vec3& T,
                        const Vec3& B,
                        const Vec3& specularColor,
                        const Vec3& diffuseColor,
                        double      shininessX,
                        double      shininessY);

    Vec3 getValue(const Vec3& inDir, const Vec3& outDir) const
    {
        const Vec3 N = Vec3(0.0, 0.0, 1.0);
        const Vec3 T = Vec3(1.0, 0.0, 0.0);
        const Vec3 B = Vec3(0.0, -1.0, 0.0);

        return compute(inDir, outDir, N, T, B,
                       specularColor_, diffuseColor_, shininessX_, shininessY_);
    }

    bool isIsotropic() const { return false; }

    std::string getName() const { return "Ashikhmin-Shirley"; }

    std::string getDescription() const
    {
        std::string reference("Michael Ashikhmin and Peter Shirley, \"An anisotropic phong BRDF model,\" Journal of Graphics Tools, Volume 5, Issue 2, pp. 25-32, December 2000.");
        return reference;
    }

private:
    Vec3   specularColor_;
    Vec3   diffuseColor_;
    double shininessX_;
    double shininessY_;
};

/*
 * Implementation
 */

inline Vec3 AshikhminShirley::compute(const Vec3& L,
                                      const Vec3& V,
                                      const Vec3& N,
                                      const Vec3& T,
                                      const Vec3& B,
                                      const Vec3& specularColor,
                                      const Vec3& diffuseColor,
                                      double      shininessX,
                                      double      shininessY)
{
    using std::max;
    using std::pow;
    using std::sqrt;

    double dotLN = L.dot(N);
    double dotVN = V.dot(N);

    Vec3 H = (L + V).normalized();
    double dotHL = H.dot(L);
    double dotHN = H.dot(N);
    double dotHT = H.dot(T);
    double dotHB = H.dot(B);

    double sqDotHT = dotHT * dotHT;
    double sqDotHB = dotHB * dotHB;
    double sqDotHN = dotHN * dotHN;

    // specular component
    Vec3 F = computeSchlickFresnel(dotHL, specularColor);
    Vec3 sBrdf =
        sqrt((shininessX + 1) * (shininessY + 1)) / (8 * PI_D) *
        pow(dotHN, (shininessX * sqDotHT + shininessY * sqDotHB) / max(1 - sqDotHN, EPSILON_D)) /
        (dotHL * max(dotLN, dotVN)) * F;

    // diffuse component
    Vec3 dCoeff = 28 * diffuseColor / (23 * PI_D);
    Vec3 sCoeff = Vec3(1.0, 1.0, 1.0) - specularColor;
    Vec3 dBrdf = dCoeff.cwiseProduct(sCoeff)
               * (1 - pow(1 - dotLN / 2, 5))
               * (1 - pow(1 - dotVN / 2, 5));

    return sBrdf + dBrdf;
}

} // namespace lb

#endif // LIBBSDF_ASHIKHMIN_SHIRLEY_H
