// =================================================================== //
// Copyright (C) 2026 Kimura Ryo                                       //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#ifndef LIBBSDF_SHEEN_H
#define LIBBSDF_SHEEN_H

#include <libbsdf/Common/Global.h>
#include <libbsdf/ReflectanceModel/ReflectanceModel.h>

namespace lb {

/*! Sheen BRDF model. */
class Sheen : public ReflectanceModel
{
public:
    Sheen(const Vec3& color, double roughness, bool terminatorSoftening = false)
        : color_(color), roughness_(roughness), terminatorSoftening_(terminatorSoftening)
    {
        parameters_.push_back(Parameter("Color", &color_));
        parameters_.push_back(Parameter("Roughness", &roughness_, 0.01, 1.0));
        parameters_.push_back(Parameter("Terminator softening", &terminatorSoftening_));
    }

    static Vec3 compute(const Vec3& L,
                        const Vec3& V,
                        const Vec3& N,
                        const Vec3& color,
                        double      roughness,
                        bool        terminatorSoftening);

    Vec3 getValue(const Vec3& inDir, const Vec3& outDir) const
    {
        const Vec3 N = Vec3(0.0, 0.0, 1.0);
        return compute(inDir, outDir, N, color_, roughness_, terminatorSoftening_);
    }

    bool isIsotropic() const { return true; }

    std::string getName() const { return "Sheen"; }

    std::string getDescription() const
    {
        std::string reference("Alejandro Conty Estevez and Christopher Kulla, \"Production Friendly Microfacet Sheen BRDF,\" SIGGRAPH 2017 Course Notes, 2017.");
        return reference;
    }

private:
    // Computes the shadowing term G.
    static double computeG(double dotNL, double dotNV, double alpha, bool terminatorSoftening);

    Vec3   color_;
    double roughness_;
    bool   terminatorSoftening_;
};

/*
 * Implementation
 */

inline Vec3 Sheen::compute(const Vec3& L,
                           const Vec3& V,
                           const Vec3& N,
                           const Vec3& color,
                           double      roughness,
                           bool        terminatorSoftening)
{
    using std::max;
    using std::pow;

    double dotNL = max(N.dot(L), EPSILON_D);
    double dotNV = max(N.dot(V), EPSILON_D);

    Vec3 H = (L + V).normalized();
    double dotHN = N.dot(H);

    double alpha = max(roughness * roughness, EPSILON_D);
    double invAlpha = 1.0 / alpha;

    // Microfacet distribution (NDF)
    double cos2h = dotHN * dotHN;
    double sin2h = max(0.0, 1.0 - cos2h);
    double D = (2.0 + invAlpha) * pow(sin2h, invAlpha * 0.5) / (2.0 * PI_D);

    double G = computeG(dotNL, dotNV, alpha, terminatorSoftening);

    return color * D * G;
}

double Sheen::computeG(double dotNL, double dotNV, double alpha, bool terminatorSoftening)
{
    using std::exp;
    using std::pow;

    double w = (1.0 - alpha) * (1.0 - alpha);
    double a = (1.0 - w) * 21.5473 + w * 25.3245;
    double b = (1.0 - w) * 3.82987 + w * 3.32435;
    double c = (1.0 - w) * 0.19823 + w * 0.16801;
    double d = (1.0 - w) * (-1.97760) + w * (-1.27393);
    double e = (1.0 - w) * (-4.32054) + w * (-4.85967);

    auto L = [&](double x) -> double { return a / (1.0 + b * pow(x, c)) + d * x + e; };

    auto lambda = [&](double cosTheta, bool softening) -> double {
        double result;
        if (cosTheta < 0.5) {
            result = exp(L(cosTheta));
        }
        else {
            result = exp(2.0 * L(0.5) - L(1.0 - cosTheta));
        }

        if (softening) {
            // The paper's artistic tweak is non-reciprocal.
            result = pow(result, 1.0 + 2.0 * pow(1.0 - cosTheta, 8.0));
        }

        return result;
    };

    double lambdaV = lambda(dotNV, false);
    double lambdaL = lambda(dotNL, terminatorSoftening);
    return 1.0 / ((1.0 + lambdaV + lambdaL) * (4.0 * dotNV * dotNL));
}

} // namespace lb

#endif // LIBBSDF_SHEEN_H