// =================================================================== //
// Copyright (C) 2021 Kimura Ryo                                       //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#ifndef LIBBSDF_UNREAL_ENGINE_4_H
#define LIBBSDF_UNREAL_ENGINE_4_H

#include <libbsdf/Common/Global.h>
#include <libbsdf/ReflectanceModel/Fresnel.h>
#include <libbsdf/ReflectanceModel/Ggx.h>
#include <libbsdf/ReflectanceModel/ReflectanceModel.h>

namespace lb {

/*! Unreal Engine 4 reflectance model. */
class UnrealEngine4 : public ReflectanceModel
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    UnrealEngine4(const Vec3&   color,
                  float         metallic,
                  float         specular,
                  float         roughness)
                  : color_      (color),
                    metallic_   (metallic),
                    specular_   (specular),
                    roughness_  (roughness)
    {
        parameters_.push_back(Parameter("Base color",   &color_));
        parameters_.push_back(Parameter("Metallic",     &metallic_, 0.0f, 1.0f));
        parameters_.push_back(Parameter("Specular",     &specular_, 0.0f, 1.0f,
                                        "Specular reflectance at normal incidence. 0.5 corresponds to 4% reflectance."));
        parameters_.push_back(Parameter("Roughness",    &roughness_, 0.01f, 1.0f));
    }

    template <typename Vec3T, typename ColorT, typename ScalarT>
    static ColorT compute(const Vec3T&      L,
                          const Vec3T&      V,
                          const Vec3T&      N,
                          const ColorT&     color,
                          const ScalarT&    metallic,
                          const ScalarT&    specular,
                          const ScalarT&    roughness);

    Vec3 getValue(const Vec3& inDir, const Vec3& outDir) const override
    {
        const Vec3 N = Vec3(0.0, 0.0, 1.0);
        return compute(inDir, outDir, N, color_, metallic_, specular_, roughness_);
    }

    bool isIsotropic() const override { return true; }

    std::string getName() const override { return "Unreal Engine 4"; }

    static std::string getReference()
    {
        return "Brian Karis, \"Real Shading in Unreal Engine 4,\" 2013.";
    }

    std::string getDescription() const override
    {
        return "Reference: " + getReference() + "\nImplementation: https://github.com/EpicGames/UnrealEngine";
    }

    /*! Modified Schlick's approximation of Fresnel reflection. */
    template <typename ScalarT, typename ColorT>
    static ColorT computeF(const ScalarT& dotVH, const ColorT& R0);

private:
    /*! Christophe Schlick, "An Inexpensive BRDF Model for Physically-Based Rendering," 1994. */
    template <typename T>
    static T computeSchlickVis(const T& dotNL, const T& dotNV, const T& alpha);

    /*!
     * Approximation of joint Smith term for GGX.
     * Eric Heitz, "Understanding the Masking-Shadowing Function in Microfacet-Based BRDFs," 2014.
     */
    template <typename T>
    static T computeSmithJointApproxVis(const T& dotNL, const T& dotNV, const T& alpha);

    Vec3    color_;
    float   metallic_;
    float   specular_;
    float   roughness_;
};

/*
 * Implementation
 */

template <typename Vec3T, typename ColorT, typename ScalarT>
ColorT UnrealEngine4::compute(const Vec3T&      L,
                              const Vec3T&      V,
                              const Vec3T&      N,
                              const ColorT&     color,
                              const ScalarT&    metallic,
                              const ScalarT&    specular,
                              const ScalarT&    roughness)
{
    using ColorScalar = typename ColorT::Scalar;

    ScalarT dotNL = static_cast<ScalarT>(N.dot(L));
    ScalarT dotNV = static_cast<ScalarT>(N.dot(V));

    Vec3T H = (L + V).normalized();

    ScalarT dotNH = static_cast<ScalarT>(N.dot(H));
    ScalarT dotVH = static_cast<ScalarT>(V.dot(H));

    ScalarT alpha = roughness * roughness;
    ScalarT sqAlpha = alpha * alpha;

    const ColorScalar s0 = ColorScalar(0.08);
    ColorT dielectricSpecularF0 = ColorT(s0, s0, s0) * specular;
    ColorT F0 = lerp(dielectricSpecularF0, color, metallic);

    ColorT F = computeF(dotVH, F0);

    ScalarT D = Ggx::computeD(dotNH, sqAlpha);
    ScalarT Vis = computeSmithJointApproxVis(dotNL, dotNV, alpha);

    ColorT fs = F * D * Vis;
    ColorT fd = color / PI_D;

    return fd * (ScalarT(1) - metallic) + fs;
}

template <typename ScalarT, typename ColorT>
ColorT UnrealEngine4::computeF(const ScalarT& dotVH, const ColorT& R0)
{
    using std::min;
    using std::pow;
    using ColorScalar = typename ColorT::Scalar;

    ColorScalar Fc = pow(ColorScalar(1) - static_cast<ColorScalar>(dotVH), ColorScalar(5));

    // If green is less than 2%, color is modified.
    ColorScalar clampedFc = min(ColorScalar(50) * R0[1], ColorScalar(1)) * Fc;

    return ColorT(clampedFc, clampedFc, clampedFc) + (ColorScalar(1) - Fc) * R0;
}

template <typename T>
T UnrealEngine4::computeSchlickVis(const T& dotNL, const T& dotNV, const T& alpha)
{
    T k = alpha * T(0.5);
    T GV = dotNV * (T(1) - k) + k;
    T GL = dotNL * (T(1) - k) + k;
    return T(0.25) / (GV * GL);
}

template <typename T>
T UnrealEngine4::computeSmithJointApproxVis(const T& dotNL, const T& dotNV, const T& alpha)
{
    T GV = dotNL * (dotNV * (T(1) - alpha) + alpha);
    T GL = dotNV * (dotNL * (T(1) - alpha) + alpha);
    return T(0.5) / (GV + GL);
}

} // namespace lb

#endif // LIBBSDF_UNREAL_ENGINE_4_H
