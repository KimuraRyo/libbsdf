// =================================================================== //
// Copyright (C) 2017 Kimura Ryo                                       //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#ifndef LIBBSDF_GGX_ANISOTROPIC_H
#define LIBBSDF_GGX_ANISOTROPIC_H

#include <libbsdf/Common/Global.h>
#include <libbsdf/ReflectanceModel/GGX.h>
#include <libbsdf/ReflectanceModel/ReflectanceModel.h>

namespace lb {

/*! GGX anisotropic reflectance model. */
class GgxAnisotropic : public ReflectanceModel
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    GgxAnisotropic(const Vec3&  specularColor,
                   const Vec3&  diffuseColor,
                   float        roughnessX,
                   float        roughnessY)
                   : specularColor_ (specularColor),
                     diffuseColor_  (diffuseColor),
                     roughnessX_    (roughnessX),
                     roughnessY_    (roughnessY)
    {
        parameters_.push_back(Parameter("Specular color",   &specularColor_));
        parameters_.push_back(Parameter("Roughness X",      &roughnessX_));
        parameters_.push_back(Parameter("Roughness Y",      &roughnessY_));
        parameters_.push_back(Parameter("Diffuse color",    &diffuseColor_));
    }

    static Vec3 compute(const Vec3& L,
                        const Vec3& V,
                        const Vec3& N,
                        const Vec3& T,
                        const Vec3& B,
                        const Vec3& specularColor,
                        const Vec3& diffuseColor,
                        float       roughnessX,
                        float       roughnessY);

    Vec3 getValue(const Vec3& inDir, const Vec3& outDir) const
    {
        const Vec3 N = Vec3(0.0, 0.0, 1.0);
        const Vec3 T = Vec3(1.0, 0.0, 0.0);
        const Vec3 B = Vec3(0.0, -1.0, 0.0);

        return compute(inDir, outDir, N, T, B,
                       specularColor_, diffuseColor_, roughnessX_, roughnessY_);
    }

    bool isIsotropic() const { return false; }

    std::string getName() const { return "GGX Anisotropic"; }

    std::string getDescription() const
    {
        std::string reference("Brent Burley, \"Physically based shading at Disney,\" part of \"Practical physically based shading in film and game production\", SIGGRAPH 2012 Course Notes, 2012.");
        return reference;
    }

private:
    Vec3    specularColor_;
    Vec3    diffuseColor_;
    float   roughnessX_;
    float   roughnessY_;
};

/*
 * Implementation
 */

inline Vec3 GgxAnisotropic::compute(const Vec3& L,
                                    const Vec3& V,
                                    const Vec3& N,
                                    const Vec3& T,
                                    const Vec3& B,
                                    const Vec3& specularColor,
                                    const Vec3& diffuseColor,
                                    float       roughnessX,
                                    float       roughnessY)
{
    using std::min;
    using std::pow;

    float alphaX = roughnessX * roughnessX;
    float alphaY = roughnessY * roughnessY;

    float roughnessXg = 0.5f + roughnessX * 0.5f;
    float roughnessYg = 0.5f + roughnessY * 0.5f;
    float alphaXg = roughnessXg * roughnessXg;
    float alphaYg = roughnessYg * roughnessYg;

    float dotLN = L.dot(N);
    float dotVN = V.dot(N);

    Vec3 H = (L + V).normalized();
    float dotHN = H.dot(N);
    float dotHT = H.dot(T);
    float dotHB = H.dot(B);
    float dotVH = min(V.dot(H), 1.0f);

    float sqAlpha = alphaX * alphaY;
    float sqAlphaG = alphaXg * alphaYg;

    // specular component
    float denominator = dotHT * dotHT / (alphaX * alphaX)
                      + dotHB * dotHB / (alphaY * alphaY)
                      + dotHN * dotHN;
    float D = 1.0f / (PI_F * sqAlpha * denominator * denominator);
    Vec3  F = schlickFresnel(dotVH, specularColor);
    float G = Ggx::computeG1(dotVN, sqAlphaG) * Ggx::computeG1(dotLN, sqAlphaG);
    Vec3 sBrdf = D * F * G / (4.0f * dotLN * dotVN);

    // diffuse component
    float Fd90 = 0.5f + 2.0f * (roughnessX + roughnessY) / 2.0f * dotVH * dotVH;
    Vec3 dBrdf = diffuseColor / PI_F
               * (1.0f + (Fd90 - 1.0f) * pow(1.0f - dotLN, 5.0f))
               * (1.0f + (Fd90 - 1.0f) * pow(1.0f - dotVN, 5.0f));

    return sBrdf + dBrdf;
}

} // namespace lb

#endif // LIBBSDF_GGX_ANISOTROPIC_H
