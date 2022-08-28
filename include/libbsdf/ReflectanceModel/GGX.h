// =================================================================== //
// Copyright (C) 2017-2021 Kimura Ryo                                  //
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
        float       roughness)
        : color_    (color),
          roughness_(roughness)
    {
        parameters_.push_back(Parameter("Color",        &color_));
        parameters_.push_back(Parameter("Roughness",    &roughness_, 0.01f, 1.0f));
    }

    Ggx(const Vec3& color,
        float       roughness,
        float       refractiveIndex,
        float       extinctionCoefficient)
        : color_                (color),
          roughness_            (roughness),
          refractiveIndex_      (refractiveIndex),
          extinctionCoefficient_(extinctionCoefficient)
    {
        parameters_.push_back(Parameter("Color",                    &color_));
        parameters_.push_back(Parameter("Roughness",                &roughness_, 0.01f, 1.0f));
        parameters_.push_back(Parameter("Refractive index",         &refractiveIndex_, 0.01f, 100.0f));
        parameters_.push_back(Parameter("Extinction coefficient",   &extinctionCoefficient_, 0.0f, 100.0f));
    }

    template <typename Vec3T, typename ColorT, typename ScalarT>
    static ColorT compute(const Vec3T&      L,
                          const Vec3T&      V,
                          const Vec3T&      N,
                          const ColorT&     color,
                          const ScalarT&    roughness);

    template <typename Vec3T, typename ColorT, typename ScalarT>
    static ColorT compute(const Vec3T&      L,
                          const Vec3T&      V,
                          const Vec3T&      N,
                          const ColorT&     color,
                          const ScalarT&    roughness,
                          const ScalarT&    refractiveIndex,
                          const ScalarT&    extinctionCoefficient);

    Vec3 getValue(const Vec3& inDir, const Vec3& outDir) const override
    {
        const Vec3 N = Vec3(0.0, 0.0, 1.0);

#if defined(LIBBSDF_USE_COLOR_INSTEAD_OF_REFRACTIVE_INDEX)
        return compute(inDir, outDir, N, color_, roughness_);
#else
        return compute(inDir, outDir, N, color_, roughness_, refractiveIndex_, extinctionCoefficient_);
#endif
    }

    bool isIsotropic() const override { return true; }

    std::string getName() const override { return "GGX (isotropic)"; }

    static std::string getReference()
    {
        return "Bruce Walter, Stephen R. Marschner, Hongsong Li, and Kenneth E. Torrance, \"Microfacet models for refraction through rough surfaces,\" Eurographics Symposium on Rendering (2007), pp. 195-206, June 2007.";
    }

    std::string getDescription() const override
    {
        return "Reference: " + getReference();
    }

    template <typename T>
    static T computeG1(const T& dotN, const T& sqAlpha);

    template <typename T>
    static T computeD(const T& dotN, const T& sqAlpha);

    template <typename T>
    static T clampDotLH(const T& dotLH);

private:
    Vec3    color_;
    float   roughness_;
    float   refractiveIndex_;
    float   extinctionCoefficient_;
};

/*
 * Implementation
 */

template <typename Vec3T, typename ColorT, typename ScalarT>
ColorT Ggx::compute(const Vec3T&    L,
                    const Vec3T&    V,
                    const Vec3T&    N,
                    const ColorT&   color,
                    const ScalarT&  roughness)
{
    using std::abs;
    using std::acos;
    using std::min;

    ScalarT dotNL = static_cast<ScalarT>(N.dot(L));
    ScalarT dotNV = static_cast<ScalarT>(N.dot(V));

    Vec3T H = (L + V).normalized();

    ScalarT dotNH = static_cast<ScalarT>(N.dot(H));
    ScalarT dotLH = static_cast<ScalarT>(clampDotLH(L.dot(H)));

    ColorT F = computeSchlickFresnel(dotLH, color);

    ScalarT alpha = roughness * roughness;
    ScalarT sqAlpha = alpha * alpha;

    ScalarT G = computeG1(dotNL, sqAlpha) * computeG1(dotNV, sqAlpha);
    ScalarT D = computeD(dotNH, sqAlpha);

    return F * G * D / (ScalarT(4) * dotNL * dotNV);
}

template <typename Vec3T, typename ColorT, typename ScalarT>
ColorT Ggx::compute(const Vec3T&    L,
                    const Vec3T&    V,
                    const Vec3T&    N,
                    const ColorT&   color,
                    const ScalarT&  roughness,
                    const ScalarT&  refractiveIndex,
                    const ScalarT&  extinctionCoefficient)
{
    using std::abs;
    using std::acos;
    using std::min;

    ScalarT dotNL = static_cast<ScalarT>(N.dot(L));
    ScalarT dotNV = static_cast<ScalarT>(N.dot(V));

    bool reflected = (dotNV >= ScalarT(0));

    // If the transmission of conductor is found, 0 is returned.
    if (!reflected && extinctionCoefficient > ScalarT(0.00001)) {
        return ColorT::Zero();
    }

    // If the refractive index of dielectric is 1, 0 is returned.
    if (refractiveIndex == ScalarT(1) && extinctionCoefficient < ScalarT(0.00001)) {
        return ColorT::Zero();
    }

    Vec3T H = reflected ? (L + V).normalized() : -(L + V * refractiveIndex).normalized();

    // Transmission from the back side of the surface.
    if (!reflected && refractiveIndex < ScalarT(1)) {
        H = -H;
    }

    ScalarT dotNH = static_cast<ScalarT>(N.dot(H));
    ScalarT dotLH = static_cast<ScalarT>(clampDotLH(L.dot(H)));
    ScalarT dotVH = reflected ? dotLH : clamp(static_cast<ScalarT>(V.dot(H)), ScalarT(-1), ScalarT(1));

    if (!reflected && (dotLH < ScalarT(0) ||
                       dotLH * dotNL < ScalarT(0) ||
                       dotVH * dotNV < ScalarT(0) ||
                       dotNH < ScalarT(0))) {
        return ColorT::Zero();
    }

    ColorT F = color * computeComplexFresnel(acos(dotLH), refractiveIndex, extinctionCoefficient);

    ScalarT alpha = roughness * roughness;
    ScalarT sqAlpha = alpha * alpha;

    ScalarT G = computeG1(dotNL, sqAlpha) * computeG1(dotNV, sqAlpha);
    ScalarT D = computeD(dotNH, sqAlpha);

    if (reflected) {
        return F * G * D / (ScalarT(4) * abs(dotNL) * abs(dotNV));
    }
    else {
        ScalarT d = dotLH + refractiveIndex * dotVH;
        return (abs(dotLH) * abs(dotVH)) / (abs(dotNL) * abs(dotNV)) *
               refractiveIndex * refractiveIndex *
               (ColorT::Ones() - F) * G * D / (d * d);
    }
}

template <typename T>
T Ggx::computeG1(const T& dotN, const T& sqAlpha)
{
    using std::sqrt;

    T sqTanN = T(1) / (dotN * dotN) - T(1);
    return T(2) / (T(1) + sqrt(T(1) + sqAlpha * sqTanN));
}

template <typename T>
T Ggx::computeD(const T& dotNH, const T& sqAlpha)
{
    T d = dotNH * dotNH * (sqAlpha - T(1)) + T(1);
    return sqAlpha / (T(PI_D) * d * d);
}

template <typename T>
T Ggx::clampDotLH(const T& dotLH)
{
#if defined(LIBBSDF_USE_CERES_SOLVER)
    // ceres::acos returns NaN as the derivative if the argument is 1 or -1.
    const T esp = std::numeric_limits<T>::epsilon();
    return clamp(static_cast<T>(dotLH), T(-1) + esp, T(1) - esp);
#else
    return clamp(static_cast<T>(dotLH), T(-1), T(1));
#endif
}

} // namespace lb

#endif // LIBBSDF_GGX_H
