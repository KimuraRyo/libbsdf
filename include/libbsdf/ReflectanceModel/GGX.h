// =================================================================== //
// Copyright (C) 2017-2019 Kimura Ryo                                  //
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

/*! GGX (Trowbridge-Reitz) BSDF model. */
class Ggx : public ReflectanceModel
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    Ggx(const Vec3& color,
        float       roughness,
        float       refractiveIndex = 1.5f,
        float       extinctionCoefficient = 0.0f)
        : color_                (color),
          roughness_            (roughness),
          refractiveIndex_      (refractiveIndex),
          extinctionCoefficient_(extinctionCoefficient)
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
                        float       refractiveIndex = 1.5f,
                        float       extinctionCoefficient = 0.0f);

    Vec3 getValue(const Vec3& inDir, const Vec3& outDir) const
    {
        const Vec3 N = Vec3(0.0, 0.0, 1.0);
        return compute(inDir, outDir, N, color_, roughness_, refractiveIndex_, extinctionCoefficient_);
    }

    bool isIsotropic() const { return true; }

    std::string getName() const { return "GGX (isotropic)"; }

    static std::string getReference()
    {
        return "Bruce Walter, Stephen R. Marschner, Hongsong Li, and Kenneth E. Torrance, \"Microfacet models for refraction through rough surfaces,\" Eurographics Symposium on Rendering (2007), pp. 195-206, June 2007.";
    }

    std::string getDescription() const
    {
        return "Reference: " + getReference();
    }

    static double computeG1(double dotN, double sqAlpha);

private:
    Vec3    color_;
    float   roughness_;
    float   refractiveIndex_;
    float   extinctionCoefficient_;
};

/*
 * Implementation
 */

inline Vec3 Ggx::compute(const Vec3&    L,
                         const Vec3&    V,
                         const Vec3&    N,
                         const Vec3&    color,
                         float          roughness,
                         float          refractiveIndex,
                         float          extinctionCoefficient)
{
    using std::abs;
    using std::acos;
    using std::min;

    double dotLN = L.dot(N);
    double dotVN = V.dot(N);

#if defined(LIBBSDF_USE_COLOR_INSTEAD_OF_REFRACTIVE_INDEX)
    Vec3 H = (L + V).normalized();

    double dotHN = H.dot(N);
    double dotLH = min(L.dot(H), Vec3::Scalar(1));
    double dotVH = min(V.dot(H), Vec3::Scalar(1));

    Vec3 F = fresnelSchlick(dotVH, color);
#else
    bool reflected = (dotVN >= 0.0f);

    // If the transmission of conductor is found, 0.0 is returned.
    if (!reflected && extinctionCoefficient > 0.00001f) {
        return Vec3::Zero();
    }

    // If the refractive index of dielectric is 1.0, 0.0 is returned.
    if (refractiveIndex == 1.0f && extinctionCoefficient < 0.00001f) {
        return Vec3::Zero();
    }

    Vec3 H = reflected ? (L + V).normalized() : -(L + refractiveIndex * V).normalized();

    // incoming direction of transmission at inside of surface
    if (!reflected && refractiveIndex < 1.0f) {
        H = -H;
    }

    double dotHN = H.dot(N);
    double dotLH = clamp(static_cast<double>(L.dot(H)), -1.0, 1.0);
    double dotVH = clamp(static_cast<double>(V.dot(H)), -1.0, 1.0);

    if (!reflected && (dotLH < 0.0 || // F
                       dotLH * dotLN < 0.0 || // G
                       dotVH * dotVN < 0.0 || // G
                       dotHN < 0.0 // D
                       )) {
        return Vec3::Zero();
    }

    float inTheta = static_cast<float>(acos(dotLH));
    Vec3 F = color * fresnelComplex(inTheta, refractiveIndex, extinctionCoefficient);
#endif

    double alpha = roughness * roughness;
    double sqAlpha = alpha * alpha;

    double G = computeG1(dotLN, sqAlpha) * computeG1(dotVN, sqAlpha);

    double sqDotHN = dotHN * dotHN;
    double tanHN = sqDotHN * (sqAlpha - 1.0) + 1.0;
    double D = sqAlpha / (PI_D * (tanHN * tanHN));

#if defined(LIBBSDF_USE_COLOR_INSTEAD_OF_REFRACTIVE_INDEX)
    return F * G * D / (4.0f * dotLN * dotVN);
#else
    if (reflected) {
        return F * G * D / (4.0 * abs(dotLN) * abs(dotVN));
    }
    else {
        double denominator = dotLH + refractiveIndex * dotVH;
        return (abs(dotLH) * abs(dotVH)) / (abs(dotLN) * abs(dotVN)) *
               refractiveIndex * refractiveIndex *
               (Vec3::Ones() - F) * G * D / (denominator * denominator);
    }
#endif
}

inline double Ggx::computeG1(double dotN, double sqAlpha)
{
    assert(sqAlpha > 0.0);

    using std::sqrt;

    double sqTanN = 1.0 / (dotN * dotN) - 1.0;
    return 2.0 / (1.0 + sqrt(1.0 + sqAlpha * sqTanN));
}

} // namespace lb

#endif // LIBBSDF_GGX_H
