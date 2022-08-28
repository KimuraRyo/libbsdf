// =================================================================== //
// Copyright (C) 2015-2021 Kimura Ryo                                  //
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

    template <typename T>
    static T compute(const T& color);

    Vec3 getValue(const Vec3& inDir, const Vec3&) const override
    {
        const Vec3 N = Vec3(0.0, 0.0, 1.0);
        return color_ * inDir.dot(N);
    }

    Vec3 getBrdfValue(const Vec3&, const Vec3&) const override
    {
        return compute(color_);
    }

    bool isIsotropic() const override { return true; }

    std::string getName() const override { return "Lambertian"; }

private:
    Vec3 color_;
};

/*
 * Implementation
 */

template <typename T>
T Lambertian::compute(const T& color)
{
    using Scalar = typename T::Scalar;
    return color / Scalar(PI_D);
}

} // namespace lb

#endif // LIBBSDF_LAMBERTIAN_H
