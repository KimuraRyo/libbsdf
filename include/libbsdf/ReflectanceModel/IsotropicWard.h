// =================================================================== //
// Copyright (C) 2015-2023 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#ifndef LIBBSDF_ISOTROPIC_WARD_H
#define LIBBSDF_ISOTROPIC_WARD_H

#include <libbsdf/Common/Global.h>
#include <libbsdf/ReflectanceModel/ReflectanceModel.h>

namespace lb {

/*! Ward isotropic reflectance model. */
class IsotropicWard : public ReflectanceModel
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    IsotropicWard(const Vec3& color, double roughness) : color_(color), roughness_(roughness)
    {
        parameters_.push_back(Parameter("Color", &color_));
        parameters_.push_back(Parameter("Roughness", &roughness_, 0.01, 1.0));
    }

    static Vec3
    compute(const Vec3& L, const Vec3& V, const Vec3& N, const Vec3& color, double roughness);

    Vec3 getValue(const Vec3& inDir, const Vec3& outDir) const
    {
        const Vec3 N = Vec3(0.0, 0.0, 1.0);
        return compute(inDir, outDir, N, color_, roughness_);
    }

    bool isIsotropic() const { return true; }

    std::string getName() const { return "Ward (isotropic)"; }

    std::string getDescription() const
    {
        std::string reference("Gregory J. Ward, \"Measuring and modeling anisotropic reflection,\" Computer Graphics (SIGGRAPH '92 Proceedings), pp. 265-272, July 1992.");
        return reference;
    }

private:
    Vec3   color_;
    double roughness_;
};

/*
 * Implementation
 */

inline Vec3 IsotropicWard::compute(const Vec3& L,
                                   const Vec3& V,
                                   const Vec3& N,
                                   const Vec3& color,
                                   double      roughness)
{
    using std::acos;
    using std::exp;
    using std::max;
    using std::sqrt;
    using std::tan;

    double dotLN = L.dot(N);
    double dotVN = V.dot(N);

    Vec3 H = (L + V).normalized();
    double dotHN = H.dot(N);

    double sqRoughness = roughness * roughness;
    double tanHN = tan(acos(dotHN));

    constexpr double suppressionCoeff = 0.01;

    double brdf = 1 / sqrt(max(dotLN * dotVN, suppressionCoeff)) *
                  exp(-(tanHN * tanHN / sqRoughness)) / (4 * PI_D * sqRoughness);
    return color * brdf;
}

} // namespace lb

#endif // LIBBSDF_ISOTROPIC_WARD_H
