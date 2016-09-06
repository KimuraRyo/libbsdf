// =================================================================== //
// Copyright (C) 2016 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#ifndef LIBBSDF_MODIFIED_PHONG_H
#define LIBBSDF_MODIFIED_PHONG_H

#include <libbsdf/Common/Global.h>
#include <libbsdf/Common/Utility.h>
#include <libbsdf/ReflectanceModel/ReflectanceModel.h>

namespace lb {

/*! Modified Phong reflectance model. */
struct ModifiedPhong : public ReflectanceModel
{
    explicit ModifiedPhong(float shininess) : shininess_(shininess)
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

    bool isIsotropic() const { return true; }

    std::string getName() const { return "Modified phong"; }

    std::string getDescription() const
    {
        std::string reference("Eric P. Lafortune and Yves D. Willems, \"Using the modified phong reflectance model for physically based rendering,\" Technical Report CW197, Leuven, Belgium, 1994.");
        return reference;
    }

private:
    float shininess_;
};

/*
 * Implementation
 */

inline float ModifiedPhong::compute(const Vec3& L,
                                    const Vec3& V,
                                    const Vec3& N,
                                    float       shininess)
{
    using std::max;
    using std::pow;

    Vec3 R = reflect(L, N);
    float dotRV = R.dot(V);
    return (shininess + 2.0f) / (2.0f * PI_F) * pow(max(dotRV, 0.0f), shininess);
}

} // namespace lb

#endif // LIBBSDF_MODIFIED_PHONG_H
