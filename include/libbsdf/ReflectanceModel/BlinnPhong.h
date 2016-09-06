// =================================================================== //
// Copyright (C) 2016 Kimura Ryo                                       //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#ifndef LIBBSDF_BLINN_PHONG_H
#define LIBBSDF_BLINN_PHONG_H

#include <libbsdf/Common/Global.h>
#include <libbsdf/ReflectanceModel/ReflectanceModel.h>

namespace lb {

/*! Blinn-Phong reflectance model. */
struct BlinnPhong : public ReflectanceModel
{
    explicit BlinnPhong(float shininess) : shininess_(shininess)
    {
        parameters_["Shininess"] = &shininess_;
    }

    static float compute(const Vec3&    L,
                         const Vec3&    V,
                         const Vec3&    N,
                         float          shininess);

    float getValue(const Vec3& inDir, const Vec3& outDir) const
    {
        const Vec3 N = Vec3(0.0, 0.0, 1.0);
        return compute(inDir, outDir, N, shininess_);
    }

    float getBrdfValue(const Vec3& inDir, const Vec3& outDir) const
    {
        using std::max;

        const Vec3 N = Vec3(0.0, 0.0, 1.0);
        float dotLN = inDir.dot(N);
        return getValue(inDir, outDir) / max(dotLN, EPSILON_F);
    }

    bool isIsotropic() const { return true; }

    std::string getName() const { return "Blinn-Phong"; }

    std::string getDescription() const
    {
        std::string reference("James F. Blinn, \"Models of light reflection for computer synthesized pictures,\" ACM Computer Graphics (SIGGRAPH '77 Proceedings), pp. 192-198, July 1977.");
        return reference;
    }

private:
    float shininess_;
};

/*
 * Implementation
 */

inline float BlinnPhong::compute(const Vec3&    L,
                                 const Vec3&    V,
                                 const Vec3&    N,
                                 float          shininess)
{
    using std::max;
    using std::pow;

    Vec3 H = (L + V).normalized();
    float dotHN = H.dot(N);
    return pow(max(dotHN, 0.0f), shininess);
}

} // namespace lb

#endif // LIBBSDF_BLINN_PHONG_H
