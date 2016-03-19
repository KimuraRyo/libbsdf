// =================================================================== //
// Copyright (C) 2015-2016 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#ifndef LIBBSDF_WARD_ISOTROPIC_H
#define LIBBSDF_WARD_ISOTROPIC_H

#include <libbsdf/Common/Global.h>
#include <libbsdf/ReflectanceModel/ReflectanceModel.h>

namespace lb {

/*! Ward isotropic reflectance model. */
struct WardIsotropic : public ReflectanceModel
{
    explicit WardIsotropic(float roughness) : roughness_(roughness)
    {
        parameters_["Roughness"] = &roughness_;
    }

    static float getResult(const Vec3&  inDir,
                           const Vec3&  outDir,
                           const Vec3&  normalDir,
                           float        roughness);

    float getValue(const Vec3& inDir, const Vec3& outDir) const
    {
        const Vec3 N = Vec3(0.0, 0.0, 1.0);
        return getResult(inDir, outDir, N, roughness_);
    }

    bool isIsotropic() const { return true; }

    std::string getName() const { return "Ward isotropic"; }

private:
    float roughness_;
};

/*
 * Implementation
 */

inline float WardIsotropic::getResult(const Vec3&   inDir,
                                      const Vec3&   outDir,
                                      const Vec3&   normalDir,
                                      float         roughness)
{
    using std::acos;
    using std::exp;
    using std::sqrt;
    using std::tan;

    float dotLN = inDir.dot(normalDir);
    float dotVN = outDir.dot(normalDir);

    Vec3 H = (inDir + outDir).normalized();
    float dotHN = H.dot(normalDir);

    float sqRoughness = roughness * roughness;
    float tanHN = tan(acos(dotHN));

    float brdf = 1.0f / sqrt(dotLN * dotVN)
               * exp(-(tanHN * tanHN / sqRoughness))
               / (4.0f * PI_F * sqRoughness);
    return brdf;
}

} // namespace lb

#endif // LIBBSDF_WARD_ISOTROPIC_H
