// =================================================================== //
// Copyright (C) 2017 Kimura Ryo                                       //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#ifndef LIBBSDF_GGX_H
#define LIBBSDF_GGX_H

#include <libbsdf/Common/Global.h>
#include <libbsdf/ReflectanceModel/Fresnel.h>
#include <libbsdf/ReflectanceModel/ReflectanceModel.h>

namespace lb {

/*! GGX (Trowbridge-Reitz) reflectance model. */
class Ggx : public ReflectanceModel
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    Ggx(const Vec3& color,
        float       roughness)
        : color_    (color),
          roughness_(roughness)
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

    std::string getName() const { return "GGX (Trowbridge-Reitz)"; }

    std::string getDescription() const
    {
        std::string reference("Bruce Walter, Stephen R. Marschner, Hongsong Li, and Kenneth E. Torrance, \"Microfacet models for refraction through rough surfaces,\" Eurographics Symposium on Rendering (2007), pp. 195-206, June 2007.");
        return reference;
    }

    static float computeG1(float dotN, float sqAlpha);

private:
    Vec3    color_;
    float   roughness_;
};

/*
 * Implementation
 */

inline Vec3 Ggx::compute(const Vec3&    L,
                         const Vec3&    V,
                         const Vec3&    N,
                         const Vec3&    color,
                         float          roughness)
{
    using std::min;

    float alpha = roughness * roughness;

    float dotLN = L.dot(N);
    float dotVN = V.dot(N);

    Vec3 H = (L + V).normalized();
    float dotHN = H.dot(N);
    float dotVH = min(V.dot(H), 1.0f);

    float sqDotHN = dotHN * dotHN;
    float sqAlpha = alpha * alpha;
    float tanHN = sqDotHN * (sqAlpha - 1.0f) + 1.0f;

    float D = sqAlpha / (PI_F * (tanHN * tanHN));
    Vec3  F = schlickFresnel(dotVH, color);
    float G = computeG1(dotVN, sqAlpha) * computeG1(dotLN, sqAlpha);

    return D * F * G / (4.0f * dotLN * dotVN);
}

inline float Ggx::computeG1(float dotN, float sqAlpha)
{
    using std::sqrt;

    return 2.0f * dotN / (dotN + sqrt(sqAlpha + (1.0f - sqAlpha) * dotN * dotN));
}

} // namespace lb

#endif // LIBBSDF_GGX_H
