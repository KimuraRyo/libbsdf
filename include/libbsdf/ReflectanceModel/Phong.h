// =================================================================== //
// Copyright (C) 2015-2019 Kimura Ryo                                  //
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
class Phong : public ReflectanceModel
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    Phong(const Vec3&   color,
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

    Vec3 getBrdfValue(const Vec3& inDir, const Vec3& outDir) const
    {
        using std::max;

        const Vec3 N = Vec3(0.0, 0.0, 1.0);
        float dotLN = inDir.dot(N);
        return getValue(inDir, outDir) / max(dotLN, EPSILON_F);
    }

    bool isIsotropic() const { return true; }

    std::string getName() const { return "Phong"; }

    std::string getDescription() const
    {
        std::string reference("Bui Tuong Phong, \"Illumination for computer generated pictures,\" Communications of the ACM, vol. 18, no. 6, pp. 311-317, June 1975.");
        return reference;
    }

private:
    Vec3    color_;
    float   shininess_;
};

/*
 * Implementation
 */

inline Vec3 Phong::compute(const Vec3&  L,
                           const Vec3&  V,
                           const Vec3&  N,
                           const Vec3&  color,
                           float        shininess)
{
    using std::max;
    using std::pow;

    Vec3 R = reflect(L, N);
    return color * pow(max(R.dot(V), Vec3::Scalar(0)), shininess);
}

} // namespace lb

#endif // LIBBSDF_PHONG_H
