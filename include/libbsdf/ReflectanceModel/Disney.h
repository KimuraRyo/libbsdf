// =================================================================== //
// Copyright (C) 2018-2019 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#ifndef LIBBSDF_DISNEY_H
#define LIBBSDF_DISNEY_H

#include <libbsdf/Common/Global.h>
#include <libbsdf/ReflectanceModel/GGX.h>
#include <libbsdf/ReflectanceModel/ReflectanceModel.h>

namespace lb {

/*! Disney reflectance model by Brent Burley et al. */
class Disney : public ReflectanceModel
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        Disney(const Vec3&  specularColor,
               const Vec3&  diffuseColor,
               float        roughnessX,
               float        roughnessY)
               : specularColor_ (specularColor),
                 diffuseColor_  (diffuseColor),
                 roughnessX_    (roughnessX),
                 roughnessY_    (roughnessY)
    {
        parameters_.push_back(Parameter("Specular color",   &specularColor_));
        parameters_.push_back(Parameter("Roughness X",      &roughnessX_, 0.01f, 1.0f));
        parameters_.push_back(Parameter("Roughness Y",      &roughnessY_, 0.01f, 1.0f));
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

    std::string getName() const { return "Disney"; }

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

inline Vec3 Disney::compute(const Vec3& L,
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

    double alphaX = roughnessX * roughnessX;
    double alphaY = roughnessY * roughnessY;
    double sqAlpha = alphaX * alphaY;

    double dotLN = L.dot(N);
    double dotVN = V.dot(N);

    Vec3 H = (L + V).normalized();

    double dotHN = H.dot(N);
    double dotHT = H.dot(T);
    double dotHB = H.dot(B);
    double dotVH = min(V.dot(H), Vec3::Scalar(1));

    Vec3 F = fresnelSchlick(dotVH, specularColor);

    // Remap roughness.
    double roughnessXG = 0.5 + roughnessX * 0.5;
    double roughnessYG = 0.5 + roughnessY * 0.5;
    double alphaXG = roughnessXG * roughnessXG;
    double alphaYG = roughnessYG * roughnessYG;
    double sqAlphaG = alphaXG * alphaYG;
    double G = Ggx::computeG1(dotVN, sqAlphaG) * Ggx::computeG1(dotLN, sqAlphaG);

    double denominatorD = dotHT * dotHT / (alphaX * alphaX)
                        + dotHB * dotHB / (alphaY * alphaY)
                        + dotHN * dotHN;
    double D = 1.0 / (PI_D * sqAlpha * denominatorD * denominatorD);

    // specular component
    Vec3 sBrdf = F * G * D / (4.0 * dotLN * dotVN);

    double Fd90 = 0.5 + 2.0 * (roughnessX + roughnessY) / 2.0 * dotVH * dotVH;

    // diffuse component
    Vec3 dBrdf = diffuseColor / PI_D
               * (1.0 + (Fd90 - 1.0) * pow(1.0 - dotLN, 5.0))
               * (1.0 + (Fd90 - 1.0) * pow(1.0 - dotVN, 5.0));

    return sBrdf + dBrdf;
}

} // namespace lb

#endif // LIBBSDF_DISNEY_H
