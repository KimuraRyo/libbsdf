// =================================================================== //
// Copyright (C) 2016-2018 Kimura Ryo                                  //
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
class ModifiedPhong : public ReflectanceModel
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    ModifiedPhong(const Vec3&   color,
                  float         shininess)
                  : color_      (color),
                    shininess_  (shininess)
    {
        parameters_.push_back(Parameter("Color",        &color_));
        parameters_.push_back(Parameter("Shininess",    &shininess_, 0.0f, 1000.0f));
    }

    static Vec3 compute(const Vec3& L,
                        const Vec3& V,
                        const Vec3& N,
                        const Vec3& color,
                        float       shininess);

    Vec3 getValue(const Vec3& inDir, const Vec3& outDir) const
    {
        const Vec3 N = Vec3(0.0, 0.0, 1.0);
        return compute(inDir, outDir, N, color_, shininess_);
    }

    bool isIsotropic() const { return true; }

    std::string getName() const { return "Modified Phong"; }

    std::string getDescription() const
    {
        std::string reference("Eric P. Lafortune and Yves D. Willems, \"Using the modified phong reflectance model for physically based rendering,\" Technical Report CW197, Leuven, Belgium, 1994.");
        return reference;
    }

private:
    Vec3    color_;
    float   shininess_;
};

/*
 * Implementation
 */

inline Vec3 ModifiedPhong::compute(const Vec3&  L,
                                   const Vec3&  V,
                                   const Vec3&  N,
                                   const Vec3&  color,
                                   float        shininess)
{
    using std::max;
    using std::pow;

    Vec3 R = reflect(L, N);
    float dotRV = R.dot(V);
    return color * (shininess + 2.0f) / (2.0f * PI_F) * pow(max(dotRV, 0.0f), shininess);
}

} // namespace lb

#endif // LIBBSDF_MODIFIED_PHONG_H
