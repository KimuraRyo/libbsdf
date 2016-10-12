// =================================================================== //
// Copyright (C) 2015-2016 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#ifndef LIBBSDF_COOK_TORRANCE_H
#define LIBBSDF_COOK_TORRANCE_H

#include <libbsdf/Common/Global.h>
#include <libbsdf/ReflectanceModel/Fresnel.h>
#include <libbsdf/ReflectanceModel/ReflectanceModel.h>

namespace lb {

/*! Cook-Torrance reflectance model. */
class CookTorrance : public ReflectanceModel
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    CookTorrance(const Vec3&    color,
                 float          roughness,
                 float          refractiveIndex)
                 : color_(color),
                   roughness_(roughness),
                   refractiveIndex_(refractiveIndex)
    {
        parameters_.push_back(Parameter("Color",            &color_));
        parameters_.push_back(Parameter("Roughness",        &roughness_));
        parameters_.push_back(Parameter("Refractive index", &refractiveIndex_));
    }

    static Vec3 compute(const Vec3& L,
                        const Vec3& V,
                        const Vec3& N,
                        const Vec3& color,
                        float       roughness,
                        float       refractiveIndex);
    
    Vec3 getValue(const Vec3& inDir, const Vec3& outDir) const
    {
        const Vec3 N = Vec3(0.0, 0.0, 1.0);
        return compute(inDir, outDir, N, color_, roughness_, refractiveIndex_);
    }

    bool isIsotropic() const { return true; }

    std::string getName() const { return "Cook-Torrance"; }

    std::string getDescription() const
    {
        std::string reference("Robert L. Cook and Kenneth E. Torrance, \"A reflectance model for computer graphics,\" Computer Graphics (SIGGRAPH '81 Proceedings), pp. 307-316, July 1981.");
        return reference;
    }

private:
    Vec3    color_;
    float   roughness_;
    float   refractiveIndex_;
};

/*
 * Implementation
 */

inline Vec3 CookTorrance::compute(const Vec3&   L,
                                  const Vec3&   V,
                                  const Vec3&   N,
                                  const Vec3&   color,
                                  float         roughness,
                                  float         refractiveIndex)
{
    using std::acos;
    using std::exp;
    using std::min;

    float dotLN = L.dot(N);
    float dotVN = V.dot(N);

    Vec3 H = (L + V).normalized();
    float dotHN = H.dot(N);

    float dotVH = min(V.dot(H), 1.0f);

    float sqDotHN = dotHN * dotHN;
    float sqRoughness = roughness * roughness;
    float sqTanHN = (1.0f - sqDotHN) / (sqRoughness * sqDotHN);

    float D = exp(-sqTanHN) / (4.0f * sqRoughness * sqDotHN * sqDotHN);
    float F = fresnelReflection(acos(dotVH), refractiveIndex);
    float G = min(dotHN * dotVN / dotVH,
                  dotHN * dotLN / dotVH);
    G = min(1.0f, 2.0f * G);

    return color * (1.0f / PI_F) * D * F * G / (dotLN * dotVN);
}

} // namespace lb

#endif // LIBBSDF_COOK_TORRANCE_H
