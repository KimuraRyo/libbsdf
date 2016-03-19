// =================================================================== //
// Copyright (C) 2015-2016 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#ifndef LIBBSDF_PHONG_H
#define LIBBSDF_PHONG_H

#include <libbsdf/Common/Global.h>
#include <libbsdf/Common/Utility.h>
#include <libbsdf/ReflectanceModel/ReflectanceModel.h>

namespace lb {

/*! Phong reflectance model. */
struct Phong : public ReflectanceModel
{
    explicit Phong(float shininess) : shininess_(shininess)
    {
        parameters_["Shininess"] = &shininess_;
    }

    static float getResult(const Vec3&  inDir,
                           const Vec3&  outDir,
                           const Vec3&  normalDir,
                           float        shininess);

    float getValue(const Vec3& inDir, const Vec3& outDir) const
    {
        const Vec3 N = Vec3(0.0, 0.0, 1.0);
        return getResult(inDir, outDir, N, shininess_);
    }

    float getBrdfValue(const Vec3& inDir, const Vec3& outDir) const
    {
        using std::max;

        const Vec3 N = Vec3(0.0, 0.0, 1.0);
        float dotLN = inDir.dot(N);
        return getValue(inDir, outDir) / max(dotLN, EPSILON_F);
    }

    bool isIsotropic() const { return true; }

    std::string getName() const { return "Phong"; }

private:
    float shininess_;
};

/*
 * Implementation
 */

inline float Phong::getResult(const Vec3&   inDir,
                              const Vec3&   outDir,
                              const Vec3&   normalDir,
                              float         shininess)
{
    using std::max;
    using std::pow;

    Vec3 R = reflect(inDir, normalDir);
    return pow(max(R.dot(outDir), 0.0f), shininess);
}

} // namespace lb

#endif // LIBBSDF_PHONG_H
