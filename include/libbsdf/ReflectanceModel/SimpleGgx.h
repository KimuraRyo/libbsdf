// =================================================================== //
// Copyright (C) 2022-2023 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#ifndef LIBBSDF_SIMPLE_GGX_H
#define LIBBSDF_SIMPLE_GGX_H

#include <libbsdf/Common/Global.h>
#include <libbsdf/ReflectanceModel/Fresnel.h>
#include <libbsdf/ReflectanceModel/Ggx.h>
#include <libbsdf/ReflectanceModel/ReflectanceModel.h>

namespace lb {

/*! GGX (Trowbridge-Reitz) BRDF model without the parameters of refractive index. */
class SimpleGgx : public ReflectanceModel
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    SimpleGgx(const Vec3& color, double roughness) : color_(color), roughness_(roughness)
    {
        parameters_.push_back(Parameter("Color", &color_));
        parameters_.push_back(Parameter("Roughness", &roughness_, 0.01, 1.0));
    }

    template <typename Vec3T, typename ColorT, typename ScalarT>
    static ColorT compute(const Vec3T&   L,
                          const Vec3T&   V,
                          const Vec3T&   N,
                          const ColorT&  color,
                          const ScalarT& roughness);

    Vec3 getValue(const Vec3& inDir, const Vec3& outDir) const override
    {
        const Vec3 N = Vec3(0.0, 0.0, 1.0);
        return compute(inDir, outDir, N, color_, roughness_);
    }

    bool isIsotropic() const override { return true; }

    std::string getName() const override { return "GGX (isotropic, no refractive index)"; }

    static std::string getReference()
    {
        return "Bruce Walter, Stephen R. Marschner, Hongsong Li, and Kenneth E. Torrance, \"Microfacet models for refraction through rough surfaces,\" Eurographics Symposium on Rendering (2007), pp. 195-206, June 2007.";
    }

    std::string getDescription() const override
    {
        return "Reference: " + getReference();
    }

private:
    Vec3   color_;
    double roughness_;
};

/*
 * Implementation
 */

template <typename Vec3T, typename ColorT, typename ScalarT>
ColorT SimpleGgx::compute(const Vec3T&   L,
                          const Vec3T&   V,
                          const Vec3T&   N,
                          const ColorT&  color,
                          const ScalarT& roughness)
{
    using std::abs;
    using std::acos;
    using std::min;

    ScalarT dotNL = static_cast<ScalarT>(N.dot(L));
    ScalarT dotNV = static_cast<ScalarT>(N.dot(V));

    Vec3T H = (L + V).normalized();

    ScalarT dotNH = static_cast<ScalarT>(N.dot(H));
    ScalarT dotLH = static_cast<ScalarT>(Ggx::clampDotLH(L.dot(H)));

    ColorT F = computeSchlickFresnel(dotLH, color);

    ScalarT alpha = roughness * roughness;
    ScalarT sqAlpha = alpha * alpha;

    ScalarT G = Ggx::computeG1(dotNL, sqAlpha) * Ggx::computeG1(dotNV, sqAlpha);
    ScalarT D = Ggx::computeD(dotNH, sqAlpha);

    return F * G * D / (ScalarT(4) * dotNL * dotNV);
}

} // namespace lb

#endif // LIBBSDF_SIMPLE_GGX_H
