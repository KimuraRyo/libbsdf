// =================================================================== //
// Copyright (C) 2015-2016 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#ifndef LIBBSDF_LAMBERTIAN_H
#define LIBBSDF_LAMBERTIAN_H

#include <libbsdf/Common/Global.h>
#include <libbsdf/ReflectanceModel/ReflectanceModel.h>

namespace lb {

/*! Lambertian model. */
class Lambertian : public ReflectanceModel
{
public:
    static float compute(const Vec3& L, const Vec3& N);

    float getValue(const Vec3& inDir, const Vec3& outDir) const
    {
        const Vec3 N = Vec3(0.0, 0.0, 1.0);
        return compute(inDir, N);
    }

    float getBrdfValue(const Vec3& inDir, const Vec3& outDir) const
    {
        return 1.0f / PI_F;
    }

    bool isIsotropic() const { return true; }

    std::string getName() const { return "Lambertian"; }
};

/*
 * Implementation
 */

inline float Lambertian::compute(const Vec3& L, const Vec3& N)
{
    return L.dot(N);
}

} // namespace lb

#endif // LIBBSDF_LAMBERTIAN_H
