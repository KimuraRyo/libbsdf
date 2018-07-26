// =================================================================== //
// Copyright (C) 2016-2018 Kimura Ryo                                  //
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
class BlinnPhong : public ReflectanceModel
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        BlinnPhong(const Vec3&  color,
                   float        shininess)
                   : color_(color),
                     shininess_(shininess)
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

    std::string getName() const { return "Blinn-Phong"; }

    std::string getDescription() const
    {
        std::string reference("James F. Blinn, \"Models of light reflection for computer synthesized pictures,\" ACM Computer Graphics (SIGGRAPH '77 Proceedings), pp. 192-198, July 1977.");
        return reference;
    }

private:
    Vec3    color_;
    float   shininess_;
};

/*
 * Implementation
 */

inline Vec3 BlinnPhong::compute(const Vec3& L,
                                const Vec3& V,
                                const Vec3& N,
                                const Vec3& color,
                                float       shininess)
{
    using std::max;
    using std::pow;

    Vec3 H = (L + V).normalized();
    float dotHN = H.dot(N);
    return color * pow(max(dotHN, 0.0f), shininess);
}

} // namespace lb

#endif // LIBBSDF_BLINN_PHONG_H
