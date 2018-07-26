// =================================================================== //
// Copyright (C) 2017-2018 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#ifndef LIBBSDF_GGX_ANISOTROPIC_H
#define LIBBSDF_GGX_ANISOTROPIC_H

#include <libbsdf/Common/Global.h>
#include <libbsdf/ReflectanceModel/GGX.h>
#include <libbsdf/ReflectanceModel/ReflectanceModel.h>

namespace lb {

/*!
 * GGX anisotropic BSDF model.
 * The default refractive index and extinction coefficient are values of aluminium at 550nm.
 */
class GgxAnisotropic : public ReflectanceModel
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    GgxAnisotropic(const Vec3&  color,
                   float        roughnessX,
                   float        roughnessY,
                   float        refractiveIndex = 0.96521f,
                   float        extinctionCoefficient = 6.3995f)
                   : color_                 (color),
                     roughnessX_            (roughnessX),
                     roughnessY_            (roughnessY),
                     refractiveIndex_       (refractiveIndex),
                     extinctionCoefficient_ (extinctionCoefficient)
    {
        parameters_.push_back(Parameter("Color",                    &color_));
        parameters_.push_back(Parameter("Roughness X",              &roughnessX_, 0.01f, 1.0f));
        parameters_.push_back(Parameter("Roughness Y",              &roughnessY_, 0.01f, 1.0f));
#if !defined(LIBBSDF_USE_COLOR_INSTEAD_OF_REFRACTIVE_INDEX)
        parameters_.push_back(Parameter("Refractive index",         &refractiveIndex_, 0.01f, 100.0f));
        parameters_.push_back(Parameter("Extinction coefficient",   &extinctionCoefficient_, 0.0f, 100.0f));
#endif
    }

    static Vec3 compute(const Vec3& L,
                        const Vec3& V,
                        const Vec3& N,
                        const Vec3& T,
                        const Vec3& B,
                        const Vec3& color,
                        float       roughnessX,
                        float       roughnessY,
                        float       refractiveIndex = 0.96521f,
                        float       extinctionCoefficient = 6.3995f);

    Vec3 getValue(const Vec3& inDir, const Vec3& outDir) const
    {
        const Vec3 N = Vec3(0.0, 0.0, 1.0);
        const Vec3 T = Vec3(1.0, 0.0, 0.0);
        const Vec3 B = Vec3(0.0, -1.0, 0.0);

        return compute(inDir, outDir, N, T, B,
                       color_, roughnessX_, roughnessY_,
                       refractiveIndex_, extinctionCoefficient_);
    }

    bool isIsotropic() const { return false; }

    std::string getName() const { return "GGX (anisotropic)"; }

    std::string getDescription() const
    {
        std::string reference("Brent Burley, \"Physically based shading at Disney,\" part of \"Practical physically based shading in film and game production\", SIGGRAPH 2012 Course Notes, 2012.");
        return reference;
    }

private:
    Vec3    color_;
    float   roughnessX_;
    float   roughnessY_;
    float   refractiveIndex_;
    float   extinctionCoefficient_;
};

/*
 * Implementation
 */

inline Vec3 GgxAnisotropic::compute(const Vec3& L,
                                    const Vec3& V,
                                    const Vec3& N,
                                    const Vec3& T,
                                    const Vec3& B,
                                    const Vec3& color,
                                    float       roughnessX,
                                    float       roughnessY,
                                    float       refractiveIndex,
                                    float       extinctionCoefficient)
{
    using std::abs;
    using std::acos;
    using std::min;

    double dotLN = L.dot(N);
    double dotVN = V.dot(N);

#if defined(LIBBSDF_USE_COLOR_INSTEAD_OF_REFRACTIVE_INDEX)
    Vec3 H = (L + V).normalized();

    double dotHN = H.dot(N);
    double dotHT = H.dot(T);
    double dotHB = H.dot(B);
    double dotLH = min(L.dot(H), 1.0f);
    double dotVH = min(V.dot(H), 1.0f);

    Vec3 F = fresnelSchlick(dotVH, color);
#else
    bool reflected = (dotVN >= 0.0);

    // If the transmission of conductor is found, 0.0 is returned.
    if (!reflected && extinctionCoefficient > 0.00001f) {
        return Vec3::Zero();
    }

    // If the refractive index of dielectric is 1.0, 0.0 is returned.
    if (refractiveIndex == 1.0f && extinctionCoefficient < 0.00001f) {
        return Vec3::Zero();
    }

    Vec3 H = reflected ? (L + V).normalized() : -(L + refractiveIndex * V).normalized();

    // incoming direction of transmission at inside of surface
    if (!reflected && refractiveIndex < 1.0f) {
        H = -H;
    }

    double dotHN = H.dot(N);
    double dotHT = H.dot(T);
    double dotHB = H.dot(B);
    double dotLH = clamp(static_cast<double>(L.dot(H)), -1.0, 1.0);
    double dotVH = clamp(static_cast<double>(V.dot(H)), -1.0, 1.0);

    if (!reflected && (dotLH < 0.0 || // F
                       dotLH * dotLN < 0.0 || // G
                       dotVH * dotVN < 0.0 || // G
                       dotHN < 0.0 // D
                       )) {
        return Vec3::Zero();
    }

    Vec3 F = color * fresnelComplex(acos(dotLH), refractiveIndex, extinctionCoefficient);
#endif

    double alphaX = roughnessX * roughnessX;
    double alphaY = roughnessY * roughnessY;
    double sqAlpha = alphaX * alphaY;

    double G = Ggx::computeG1(dotLN, sqAlpha) * Ggx::computeG1(dotVN, sqAlpha);

    // GTR (Generalized-Trowbridge-Reitz) distribution function is implemented here.
    // This function with gamma = 2 is equivalent to GGX.
    double denominatorD = dotHT * dotHT / (alphaX * alphaX)
                        + dotHB * dotHB / (alphaY * alphaY)
                        + dotHN * dotHN;
    double D = 1.0 / (PI_F * sqAlpha * denominatorD * denominatorD);

#if defined(LIBBSDF_USE_COLOR_INSTEAD_OF_REFRACTIVE_INDEX)
    return F * G * D / (4.0 * dotLN * dotVN);
#else
    if (reflected) {
        return F * G * D / (4.0 * abs(dotLN) * abs(dotVN));
    }
    else {
        double denominator = dotLH + refractiveIndex * dotVH;
        return (abs(dotLH) * abs(dotVH)) / (abs(dotLN) * abs(dotVN)) *
               refractiveIndex * refractiveIndex *
               (Vec3::Ones() - F) * G * D / (denominator * denominator);
    }
#endif
}

} // namespace lb

#endif // LIBBSDF_GGX_ANISOTROPIC_H
