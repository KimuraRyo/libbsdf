// =================================================================== //
// Copyright (C) 2015-2017 Kimura Ryo                                  //
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
                 float          roughness)
                 : color_       (color),
                   roughness_   (roughness)
    {
        parameters_.push_back(Parameter("Color",        &color_));
        parameters_.push_back(Parameter("Roughness",    &roughness_));
    }

    static Vec3 compute(const Vec3& L,
                        const Vec3& V,
                        const Vec3& N,
                        const Vec3& color,
                        float       roughness);
    
    Vec3 getValue(const Vec3& inDir, const Vec3& outDir) const
    {
        const Vec3 N = Vec3(0.0, 0.0, 1.0);
        return compute(inDir, outDir, N, color_, roughness_);
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
};

/*
 * Implementation
 */

inline Vec3 CookTorrance::compute(const Vec3&   L,
                                  const Vec3&   V,
                                  const Vec3&   N,
                                  const Vec3&   color,
                                  float         roughness)
{
    using std::acos;
    using std::exp;
    using std::min;

    float alpha = roughness * roughness;

    float dotLN = L.dot(N);
    float dotVN = V.dot(N);

    Vec3 H = (L + V).normalized();
    float dotHN = H.dot(N);
    float dotVH = min(V.dot(H), 1.0f);

    float sqDotHN = dotHN * dotHN;
    float sqAlpha = alpha * alpha;
    float sqTanHN = (1.0f - sqDotHN) / (sqAlpha * sqDotHN);

    float D = exp(-sqTanHN) / (4.0f * sqAlpha * sqDotHN * sqDotHN);
    Vec3  F = schlickFresnel(dotVH, color);
    float G = min(dotHN * dotVN / dotVH,
                  dotHN * dotLN / dotVH);
    G = min(1.0f, 2.0f * G);

    return (1.0f / PI_F) * D * F * G / (dotLN * dotVN);
}

} // namespace lb

#endif // LIBBSDF_COOK_TORRANCE_H
