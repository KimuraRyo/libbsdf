// =================================================================== //
// Copyright (C) 2015-2022 Kimura Ryo                                  //
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

    IsotropicWard(const Vec3&   color,
                  float         roughness)
                  : color_      (color),
                    roughness_  (roughness)
    {
        parameters_.push_back(Parameter("Color",        &color_));
        parameters_.push_back(Parameter("Roughness",    &roughness_, 0.01f, 1.0f));
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

    std::string getName() const { return "Ward (isotropic)"; }

    std::string getDescription() const
    {
        std::string reference("Gregory J. Ward, \"Measuring and modeling anisotropic reflection,\" Computer Graphics (SIGGRAPH '92 Proceedings), pp. 265-272, July 1992.");
        return reference;
    }

private:
    Vec3    color_;
    float   roughness_;
};

/*
 * Implementation
 */

inline Vec3 IsotropicWard::compute(const Vec3&  L,
                                   const Vec3&  V,
                                   const Vec3&  N,
                                   const Vec3&  color,
                                   float        roughness)
{
    using std::acos;
    using std::exp;
    using std::max;
    using std::sqrt;
    using std::tan;

    float dotLN = L.dot(N);
    float dotVN = V.dot(N);

    Vec3 H = (L + V).normalized();
    float dotHN = H.dot(N);

    float sqRoughness = roughness * roughness;
    float tanHN = tan(acos(dotHN));

    constexpr float suppressionCoeff = 0.01f;
    float brdf = 1.0f / sqrt(max(dotLN * dotVN, suppressionCoeff))
               * exp(-(tanHN * tanHN / sqRoughness))
               / (4.0f * PI_F * sqRoughness);
    return color * brdf;
}

} // namespace lb

#endif // LIBBSDF_ISOTROPIC_WARD_H
