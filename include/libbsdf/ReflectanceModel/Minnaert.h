// =================================================================== //
// Copyright (C) 2016-2018 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#ifndef LIBBSDF_MINNAERT_H
#define LIBBSDF_MINNAERT_H

#include <libbsdf/Common/Global.h>
#include <libbsdf/Common/Utility.h>
#include <libbsdf/ReflectanceModel/ReflectanceModel.h>

namespace lb {

/*! Minnaert reflectance model. */
class Minnaert : public ReflectanceModel
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    Minnaert(const Vec3&    albedo,
             float          darkness)
             : albedo_(albedo),
               darkness_(darkness)
    {
        parameters_.push_back(Parameter("Albedo",   &albedo_));
        parameters_.push_back(Parameter("Darkness", &darkness_, 0.0f, 100.0f));
    }

    /*! \warning The unit of a returned value is the radiance factor. */
    static Vec3 compute(const Vec3& L,
                        const Vec3& V,
                        const Vec3& N,
                        const Vec3& albedo,
                        float       darkness);

    Vec3 getValue(const Vec3& inDir, const Vec3& outDir) const
    {
        const Vec3 N = Vec3(0.0, 0.0, 1.0);
        return compute(inDir, outDir, N, albedo_, darkness_);
    }

    Vec3 getBrdfValue(const Vec3& inDir, const Vec3& outDir) const
    {
        return getValue(inDir, outDir) / PI_F;
    }

    bool isIsotropic() const { return true; }

    std::string getName() const { return "Minnaert"; }

    std::string getDescription() const
    {
        std::string reference("Marcel Minnaert, \"The reciprocity principle in lunar photometry,\" Astrophysical Journal, vol. 93, no. 6, pp. 403-410, May 1941.");
        return reference;
    }

private:
    Vec3    albedo_;
    float   darkness_;
};

/*
 * Implementation
 */

inline Vec3 Minnaert::compute(const Vec3&   L,
                              const Vec3&   V,
                              const Vec3&   N,
                              const Vec3&   albedo,
                              float         darkness)
{
    using std::pow;

    float dotLN = L.dot(N);
    float dotVN = V.dot(N);

    Vec3 val = albedo * pow(dotLN * dotVN, darkness - 1.0f);

    // Normalized values are returned.
    return (darkness + 1.0f) * 0.5f * val;
}

} // namespace lb

#endif // LIBBSDF_MINNAERT_H
