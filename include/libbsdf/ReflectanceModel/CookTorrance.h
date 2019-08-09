// =================================================================== //
// Copyright (C) 2015-2019 Kimura Ryo                                  //
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

/*!
 * Cook-Torrance reflectance model.
 * The default refractive index and extinction coefficient are values of aluminium at 550nm.
 */
class CookTorrance : public ReflectanceModel
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    CookTorrance(const Vec3&    color,
                 float          roughness,
                 float          refractiveIndex = 0.96521f,
                 float          extinctionCoefficient = 6.3995f)
                 : color_                   (color),
                   roughness_               (roughness),
                   refractiveIndex_         (refractiveIndex),
                   extinctionCoefficient_   (extinctionCoefficient)
    {
        parameters_.push_back(Parameter("Color",                    &color_));
        parameters_.push_back(Parameter("Roughness",                &roughness_, 0.01f, 1.0f));
#if !defined(LIBBSDF_USE_COLOR_INSTEAD_OF_REFRACTIVE_INDEX)
        parameters_.push_back(Parameter("Refractive index",         &refractiveIndex_, 0.01f, 100.0f));
        parameters_.push_back(Parameter("Extinction coefficient",   &extinctionCoefficient_, 0.0f, 100.0f));
#endif
    }

    static Vec3 compute(const Vec3& L,
                        const Vec3& V,
                        const Vec3& N,
                        const Vec3& color,
                        float       roughness,
                        float       refractiveIndex = 0.96521f,
                        float       extinctionCoefficient = 6.3995f);

    Vec3 getValue(const Vec3& inDir, const Vec3& outDir) const
    {
        const Vec3 N = Vec3(0.0, 0.0, 1.0);
        return compute(inDir, outDir, N, color_, roughness_, refractiveIndex_, extinctionCoefficient_);
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
    float   extinctionCoefficient_;
};

/*
 * Implementation
 */

inline Vec3 CookTorrance::compute(const Vec3&   L,
                                  const Vec3&   V,
                                  const Vec3&   N,
                                  const Vec3&   color,
                                  float         roughness,
                                  float         refractiveIndex,
                                  float         extinctionCoefficient)
{
    using std::acos;
    using std::exp;
    using std::min;

    double alpha = roughness * roughness;

    double dotLN = L.dot(N);
    double dotVN = V.dot(N);

    Vec3 H = (L + V).normalized();
    double dotHN = H.dot(N);
    double dotVH = min(V.dot(H), Vec3::Scalar(1));

#if defined(LIBBSDF_USE_COLOR_INSTEAD_OF_REFRACTIVE_INDEX)
    Vec3 F = fresnelSchlick(dotVH, color);
#else
    float inTheta = static_cast<float>(std::acos(dotVH));
    Vec3 F = color * fresnelComplex(inTheta, refractiveIndex, extinctionCoefficient);
#endif

    double G = min(dotHN * dotVN / dotVH,
                   dotHN * dotLN / dotVH);
    G = min(1.0, 2.0 * G);

    double sqDotHN = dotHN * dotHN;
    double sqAlpha = alpha * alpha;
    double sqTanHN = (1.0 - sqDotHN) / (sqAlpha * sqDotHN);
    double D = exp(-sqTanHN) / (PI_D * sqAlpha * sqDotHN * sqDotHN);

    return F * G * D / (4.0 * dotLN * dotVN);
}

} // namespace lb

#endif // LIBBSDF_COOK_TORRANCE_H
