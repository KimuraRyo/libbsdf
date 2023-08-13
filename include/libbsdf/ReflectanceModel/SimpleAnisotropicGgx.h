// =================================================================== //
// Copyright (C) 2022-2023 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#ifndef LIBBSDF_SIMPLE_ANISOTROPIC_GGX_H
#define LIBBSDF_SIMPLE_ANISOTROPIC_GGX_H

#include <libbsdf/Common/Global.h>
#include <libbsdf/ReflectanceModel/Ggx.h>
#include <libbsdf/ReflectanceModel/ReflectanceModel.h>

namespace lb {

/*! Anisotropic GGX BRDF model without the parameters of refractive index. */
class SimpleAnisotropicGgx : public ReflectanceModel
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    SimpleAnisotropicGgx(const Vec3& color, double roughnessX, double roughnessY)
        : color_(color), roughnessX_(roughnessX), roughnessY_(roughnessY)
    {
        parameters_.push_back(Parameter("Color", &color_));
        parameters_.push_back(Parameter("Roughness X", &roughnessX_, 0.01, 1.0));
        parameters_.push_back(Parameter("Roughness Y", &roughnessY_, 0.01, 1.0));
    }

    template <typename Vec3T, typename ColorT, typename ScalarT>
    static ColorT compute(const Vec3T&   L,
                          const Vec3T&   V,
                          const Vec3T&   N,
                          const Vec3T&   T,
                          const Vec3T&   B,
                          const ColorT&  color,
                          const ScalarT& roughnessX,
                          const ScalarT& roughnessY);

    Vec3 getValue(const Vec3& inDir, const Vec3& outDir) const override
    {
        const Vec3 N = Vec3(0.0, 0.0, 1.0);
        const Vec3 T = Vec3(1.0, 0.0, 0.0);
        const Vec3 B = Vec3(0.0, -1.0, 0.0);
        return compute(inDir, outDir, N, T, B, color_, roughnessX_, roughnessY_);
    }

    bool isIsotropic() const override { return false; }

    std::string getName() const override { return "GGX (anisotropic, no refractive index)"; }

    std::string getDescription() const override
    {
        std::string reference("Brent Burley, \"Physically based shading at Disney,\" part of \"Practical physically based shading in film and game production\", SIGGRAPH 2012 Course Notes, 2012.");
        return reference;
    }

private:
    Vec3   color_;
    double roughnessX_;
    double roughnessY_;
};

/*
 * Implementation
 */

template <typename Vec3T, typename ColorT, typename ScalarT>
ColorT SimpleAnisotropicGgx::compute(const Vec3T&   L,
                                     const Vec3T&   V,
                                     const Vec3T&   N,
                                     const Vec3T&   T,
                                     const Vec3T&   B,
                                     const ColorT&  color,
                                     const ScalarT& roughnessX,
                                     const ScalarT& roughnessY)
{
    using std::abs;
    using std::acos;
    using std::min;

    ScalarT dotNL = static_cast<ScalarT>(N.dot(L));
    ScalarT dotNV = static_cast<ScalarT>(N.dot(V));

    Vec3T H = (L + V).normalized();

    ScalarT dotNH = static_cast<ScalarT>(N.dot(H));
    ScalarT dotTH = static_cast<ScalarT>(T.dot(H));
    ScalarT dotBH = static_cast<ScalarT>(B.dot(H));

    ScalarT dotLH = static_cast<ScalarT>(Ggx::clampDotLH(L.dot(H)));

    ColorT F = computeSchlickFresnel(dotLH, color);

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

    return F * G * D / (ScalarT(4) * dotNL * dotNV);
}

} // namespace lb

#endif // LIBBSDF_SIMPLE_ANISOTROPIC_GGX_H
