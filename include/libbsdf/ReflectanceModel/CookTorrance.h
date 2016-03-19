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
struct CookTorrance : public ReflectanceModel
{
    CookTorrance(float roughness,
                 float refractiveIndex)
                 : roughness_(roughness),
                   refractiveIndex_(refractiveIndex)
    {
        parameters_["Roughness"] = &roughness_;
        parameters_["Refractive index"] = &refractiveIndex_;
    }

    static float getResult(const Vec3&  inDir,
                           const Vec3&  outDir,
                           const Vec3&  normalDir,
                           float        roughness,
                           float        refractiveIndex);
    
    float getValue(const Vec3& inDir, const Vec3& outDir) const
    {
        const Vec3 N = Vec3(0.0, 0.0, 1.0);
        return getResult(inDir, outDir, N, roughness_, refractiveIndex_);
    }

    bool isIsotropic() const { return true; }

    std::string getName() const { return "Cook-Torrance"; }

private:
    float roughness_;
    float refractiveIndex_;
};

/*
 * Implementation
 */

inline float CookTorrance::getResult(const Vec3&    inDir,
                                     const Vec3&    outDir,
                                     const Vec3&    normalDir,
                                     float          roughness,
                                     float          refractiveIndex)
{
    using std::acos;
    using std::exp;
    using std::min;

    float dotLN = inDir.dot(normalDir);
    float dotVN = outDir.dot(normalDir);

    Vec3 H = (inDir + outDir).normalized();
    float dotHN = H.dot(normalDir);

    float dotVH = min(outDir.dot(H), 1.0f);

    float sqDotHN = dotHN * dotHN;
    float sqRoughness = roughness * roughness;
    float sqTanHN = (1.0f - sqDotHN) / (sqRoughness * sqDotHN);

    float D = exp(-sqTanHN) / (4.0f * sqRoughness * sqDotHN * sqDotHN);
    float F = fresnelReflection(acos(dotVH), refractiveIndex);
    float G = min(dotHN * dotVN / dotVH,
                  dotHN * dotLN / dotVH);
    G = min(1.0f, 2.0f * G);

    return (1.0f / PI_F) * D * F * G / (dotLN * dotVN);
}

} // namespace lb

#endif // LIBBSDF_COOK_TORRANCE_H
