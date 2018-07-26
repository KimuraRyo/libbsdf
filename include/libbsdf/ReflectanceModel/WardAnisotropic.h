// =================================================================== //
// Copyright (C) 2015-2018 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#ifndef LIBBSDF_WARD_ANISOTROPIC_H
#define LIBBSDF_WARD_ANISOTROPIC_H

#include <libbsdf/Common/Global.h>
#include <libbsdf/ReflectanceModel/ReflectanceModel.h>

namespace lb {

/*! Ward anisotropic reflectance model. */
class WardAnisotropic : public ReflectanceModel
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    WardAnisotropic(const Vec3& color,
                    float       roughnessX,
                    float       roughnessY)
                    : color_(color),
                      roughnessX_(roughnessX),
                      roughnessY_(roughnessY)
    {
        parameters_.push_back(Parameter("Color",        &color_));
        parameters_.push_back(Parameter("Roughness X",  &roughnessX_, 0.01f, 1.0f));
        parameters_.push_back(Parameter("Roughness Y",  &roughnessY_, 0.01f, 1.0f));
    }

    static Vec3 compute(const Vec3& L,
                        const Vec3& V,
                        const Vec3& N,
                        const Vec3& T,
                        const Vec3& B,
                        const Vec3& color,
                        float       roughnessX,
                        float       roughnessY);

    Vec3 getValue(const Vec3& inDir, const Vec3& outDir) const
    {
        const Vec3 N = Vec3(0.0, 0.0, 1.0);
        const Vec3 T = Vec3(1.0, 0.0, 0.0);
        const Vec3 B = Vec3(0.0, -1.0, 0.0);

        return compute(inDir, outDir, N, T, B, color_, roughnessX_, roughnessY_);
    }

    bool isIsotropic() const { return false; }

    std::string getName() const { return "Ward (anisotropic)"; }

    std::string getDescription() const
    {
        std::string reference("Gregory J. Ward, \"Measuring and modeling anisotropic reflection,\" Computer Graphics (SIGGRAPH '92 Proceedings), pp. 265-272, July 1992.");
        return reference;
    }

private:
    Vec3    color_;
    float   roughnessX_;
    float   roughnessY_;
};

/*
 * Implementation
 */

inline Vec3 WardAnisotropic::compute(const Vec3&    L,
                                     const Vec3&    V,
                                     const Vec3&    N,
                                     const Vec3&    T,
                                     const Vec3&    B,
                                     const Vec3&    color,
                                     float          roughnessX,
                                     float          roughnessY)
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
    float dotHT = H.dot(T);
    float dotHB = H.dot(B);

    float sqDotHT = (dotHT / roughnessX) * (dotHT / roughnessX);
    float sqDotHB = (dotHB / roughnessY) * (dotHB / roughnessY);

    float brdf = 1.0f / sqrt(max(dotLN * dotVN, EPSILON_F))
               * exp(-2.0f * (sqDotHT + sqDotHB) / (1.0f + dotHN))
               / (4.0f * PI_F * roughnessX * roughnessY);

    return color * brdf;
}

} // namespace lb

#endif // LIBBSDF_WARD_ANISOTROPIC_H
