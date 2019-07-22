// =================================================================== //
// Copyright (C) 2015-2019 Kimura Ryo                                  //
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
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    explicit Lambertian(const Vec3& color) : color_(color)
    {
        parameters_.push_back(Parameter("Color", &color_));
    }

    static Vec3 compute(const Vec3& L,
                        const Vec3& N,
                        const Vec3& color);

    Vec3 getValue(const Vec3& inDir, const Vec3&) const
    {
        const Vec3 N = Vec3(0.0, 0.0, 1.0);
        return compute(inDir, N, color_);
    }

    Vec3 getBrdfValue(const Vec3&, const Vec3&) const
    {
        return color_ / PI_F;
    }

    bool isIsotropic() const { return true; }

    std::string getName() const { return "Lambertian"; }

private:
    Vec3 color_;
};

/*
 * Implementation
 */

inline Vec3 Lambertian::compute(const Vec3& L,
                                const Vec3& N,
                                const Vec3& color)
{
    return color * L.dot(N);
}

} // namespace lb

#endif // LIBBSDF_LAMBERTIAN_H
