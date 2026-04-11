// =================================================================== //
// Copyright (C) 2026 Kimura Ryo                                       //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#ifndef LIBBSDF_BRADY_H
#define LIBBSDF_BRADY_H

#include <libbsdf/Common/Global.h>
#include <libbsdf/ReflectanceModel/Fresnel.h>
#include <libbsdf/ReflectanceModel/ReflectanceModel.h>

namespace lb {

/*! Brady BRDF model A. */
class Brady : public ReflectanceModel
{
public:
    Brady(const Vec3& specularColor,
          const Vec3& diffuseColor,
          const Vec3& F0,
          double      alpha,
          double      beta)
        : specularColor_(specularColor),
          diffuseColor_(diffuseColor),
          F0_(F0),
          alpha_(alpha),
          beta_(beta)
    {
        parameters_.push_back(Parameter("Specular color", &specularColor_));
        parameters_.push_back(Parameter("Diffuse color", &diffuseColor_));
        parameters_.push_back(Parameter("F0", &F0_));
        parameters_.push_back(Parameter("alpha", &alpha_, 0.0, 3.0));
        parameters_.push_back(Parameter("beta", &beta_));
    }

    static Vec3 compute(const Vec3& L,
                        const Vec3& V,
                        const Vec3& N,
                        const Vec3& specularColor,
                        const Vec3& diffuseColor,
                        const Vec3& F0,
                        double      alpha,
                        double      beta);

    Vec3 getValue(const Vec3& inDir, const Vec3& outDir) const
    {
        const Vec3 N = Vec3(0.0, 0.0, 1.0);
        return compute(inDir, outDir, N, specularColor_, diffuseColor_, F0_, alpha_, beta_);
    }

    bool isIsotropic() const { return true; }

    std::string getName() const { return "Brady (BRDF Model A)"; }

    std::string getDescription() const
    {
        std::string reference("Adam Brady, Jason Lawrence, Pieter Peers, and Westley Weimer, \"genBRDF: discovering new analytic BRDFs with genetic programming,\" SIGGRAPH 2014 Technical Papers, 2014.");
        return "Reference: " + reference;
    }

private:
    Vec3   specularColor_;
    Vec3   diffuseColor_;
    Vec3   F0_;
    double alpha_;
    double beta_;
};

/*
 * Implementation
 */

inline Vec3 Brady::compute(const Vec3& L,
                           const Vec3& V,
                           const Vec3& N,
                           const Vec3& specularColor,
                           const Vec3& diffuseColor,
                           const Vec3& F0,
                           double      alpha,
                           double      beta)
{
    using std::max;

    Vec3 H = (L + V).normalized();

    double dotHN = N.dot(H);
    double dotHL = L.dot(H);

    // Delta is angle between H and N.
    double tanDelta = std::sqrt(max(1.0 - dotHN * dotHN, 0.0)) / max(dotHN, EPSILON_D);
    double DPrime = std::exp(-std::pow(tanDelta / (beta * beta), alpha));

    Vec3 F = computeSchlickFresnel(dotHL, F0);

    Vec3 diffuse = diffuseColor / PI_D;
    Vec3 specular = specularColor.cwiseProduct(F) * (DPrime / (4.0 * max(dotHL, EPSILON_D)));

    return diffuse + specular;
}

} // namespace lb

#endif // LIBBSDF_BRADY_H