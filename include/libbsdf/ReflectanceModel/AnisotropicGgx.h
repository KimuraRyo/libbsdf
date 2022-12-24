// =================================================================== //
// Copyright (C) 2017-2022 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#ifndef LIBBSDF_ANISOTROPIC_GGX_H
#define LIBBSDF_ANISOTROPIC_GGX_H

#include <libbsdf/Common/Global.h>
#include <libbsdf/ReflectanceModel/Ggx.h>
#include <libbsdf/ReflectanceModel/ReflectanceModel.h>

namespace lb {

/*! Anisotropic GGX BSDF model. */
class AnisotropicGgx : public ReflectanceModel
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    AnisotropicGgx(const Vec3& color,
                   float       roughnessX,
                   float       roughnessY,
                   float       refractiveIndex,
                   float       extinctionCoefficient = 0.0f)
        : color_(color),
          roughnessX_(roughnessX),
          roughnessY_(roughnessY),
          refractiveIndex_(refractiveIndex),
          extinctionCoefficient_(extinctionCoefficient)
    {
        parameters_.push_back(Parameter("Color", &color_));
        parameters_.push_back(Parameter("Roughness X", &roughnessX_, 0.01f, 1.0f));
        parameters_.push_back(Parameter("Roughness Y", &roughnessY_, 0.01f, 1.0f));
        parameters_.push_back(Parameter("Refractive index", &refractiveIndex_, 0.01f, 100.0f));
        parameters_.push_back(
            Parameter("Extinction coefficient", &extinctionCoefficient_, 0.0f, 100.0f));
    }

    template <typename Vec3T, typename ColorT, typename ScalarT>
    static ColorT compute(const Vec3T&   L,
                          const Vec3T&   V,
                          const Vec3T&   N,
                          const Vec3T&   T,
                          const Vec3T&   B,
                          const ColorT&  color,
                          const ScalarT& roughnessX,
                          const ScalarT& roughnessY,
                          const ScalarT& refractiveIndex,
                          const ScalarT& extinctionCoefficient);

    Vec3 getValue(const Vec3& inDir, const Vec3& outDir) const override
    {
        const Vec3 N = Vec3(0.0, 0.0, 1.0);
        const Vec3 T = Vec3(1.0, 0.0, 0.0);
        const Vec3 B = Vec3(0.0, -1.0, 0.0);
        return compute(inDir, outDir, N, T, B, color_, roughnessX_, roughnessY_, refractiveIndex_, extinctionCoefficient_);
    }

    bool isIsotropic() const override { return false; }

    std::string getName() const override { return "GGX (anisotropic)"; }

    std::string getDescription() const override
    {
        std::string reference("Brent Burley, \"Physically based shading at Disney,\" part of \"Practical physically based shading in film and game production\", SIGGRAPH 2012 Course Notes, 2012.");
        return reference;
    }

private:
    Vec3  color_;
    float roughnessX_;
    float roughnessY_;
    float refractiveIndex_;
    float extinctionCoefficient_;
};

/*
 * Implementation
 */

template <typename Vec3T, typename ColorT, typename ScalarT>
ColorT AnisotropicGgx::compute(const Vec3T&   L,
                               const Vec3T&   V,
                               const Vec3T&   N,
                               const Vec3T&   T,
                               const Vec3T&   B,
                               const ColorT&  color,
                               const ScalarT& roughnessX,
                               const ScalarT& roughnessY,
                               const ScalarT& refractiveIndex,
                               const ScalarT& extinctionCoefficient)
{
    using std::abs;
    using std::acos;
    using std::min;

    ScalarT dotNL = static_cast<ScalarT>(N.dot(L));
    ScalarT dotNV = static_cast<ScalarT>(N.dot(V));

    bool reflected = (dotNV >= ScalarT(0));

    // If the transmission of conductor is found, 0.0 is returned.
    if (!reflected && extinctionCoefficient > ScalarT(0.00001)) {
        return ColorT::Zero();
    }

    // If the refractive index of dielectric is 1.0, 0.0 is returned.
    if (refractiveIndex == ScalarT(1) && extinctionCoefficient < ScalarT(0.00001)) {
        return ColorT::Zero();
    }

    Vec3T H = reflected ? (L + V).normalized() : -(L + refractiveIndex * V).normalized();

    // Transmission from the back side of the surface.
    if (!reflected && refractiveIndex < ScalarT(1)) {
        H = -H;
    }

    ScalarT dotNH = static_cast<ScalarT>(N.dot(H));
    ScalarT dotTH = static_cast<ScalarT>(T.dot(H));
    ScalarT dotBH = static_cast<ScalarT>(B.dot(H));
    ScalarT dotLH = static_cast<ScalarT>(Ggx::clampDotLH(L.dot(H)));
    ScalarT dotVH = reflected ? dotLH : clamp(static_cast<ScalarT>(V.dot(H)), ScalarT(-1), ScalarT(1));

    if (!reflected && (dotLH < ScalarT(0) ||
                       dotLH * dotNL < ScalarT(0) ||
                       dotVH * dotNV < ScalarT(0) ||
                       dotNH < ScalarT(0))) {
        return ColorT::Zero();
    }

    ColorT F = color * computeComplexFresnel(acos(dotLH), refractiveIndex, extinctionCoefficient);

    ScalarT alphaX = roughnessX * roughnessX;
    ScalarT alphaY = roughnessY * roughnessY;
    ScalarT sqAlpha = alphaX * alphaY;

    ScalarT G = Ggx::computeG1(dotNL, sqAlpha) * Ggx::computeG1(dotNV, sqAlpha);

    // GTR (Generalized-Trowbridge-Reitz) distribution function is implemented here.
    // This function with gamma = 2 is equivalent to GGX.
    ScalarT d = dotTH * dotTH / (alphaX * alphaX)
              + dotBH * dotBH / (alphaY * alphaY)
              + dotNH * dotNH;
    ScalarT D = ScalarT(1) / (ScalarT(PI_D) * sqAlpha * d * d);

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

} // namespace lb

#endif // LIBBSDF_ANISOTROPIC_GGX_H
