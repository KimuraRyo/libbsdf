// =================================================================== //
// Copyright (C) 2026 Kimura Ryo                                      //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#ifndef LIBBSDF_VARYING_IRIDESCENCE_MICROFACET_H
#define LIBBSDF_VARYING_IRIDESCENCE_MICROFACET_H

#include <libbsdf/Common/Global.h>
#include <libbsdf/ReflectanceModel/ReflectanceModel.h>

namespace lb {

/*! GGX BRDF model with a thin-film interference term for varying iridescence. */
class VaryingIridescenceMicrofacet : public ReflectanceModel
{
public:
    VaryingIridescenceMicrofacet(double roughness,
                                 double filmThickness,
                                 double filmIor,
                                 double baseIor,
                                 double baseExtinctionCoefficient)
        : roughness_(roughness),
          filmThickness_(filmThickness),
          filmIor_(filmIor),
          baseIor_(baseIor),
          baseExtinctionCoefficient_(baseExtinctionCoefficient)
    {
        parameters_.push_back(Parameter("Roughness", &roughness_, 0.01, 1.0));
        parameters_.push_back(Parameter("Film thickness", &filmThickness_, 0.0, 10000.0,
                                        "Thickness of the thin film in nanometers"));
        parameters_.push_back(Parameter("Film refractive index", &filmIor_, 1.0, 5.0));
        parameters_.push_back(Parameter("Base refractive index", &baseIor_, 1.0, 5.0));
        parameters_.push_back(Parameter("Base extinction coefficient", &baseExtinctionCoefficient_,
                                        0.0, 5.0,
                                        "0 for a dielectric base, > 0 for a conductor base"));
    }

    /*! Evaluates a BRDF value with the interference of light reflected by the thin film. */
    static Vec3 compute(const Vec3& L,
                        const Vec3& V,
                        const Vec3& N,
                        double      roughness,
                        double      filmThickness,
                        double      filmIor,
                        double      baseIor,
                        double      baseExtinctionCoefficient);

    Vec3 getValue(const Vec3& inDir, const Vec3& outDir) const override
    {
        const Vec3 N = Vec3(0.0, 0.0, 1.0);
        return compute(inDir, outDir, N, roughness_, filmThickness_, filmIor_, baseIor_,
                       baseExtinctionCoefficient_);
    }

    bool isIsotropic() const override { return true; }

    std::string getName() const override { return "Varying iridescence microfacet"; }

    static std::string getReference()
    {
        return "Laurent Belcour and Pascal Barla, \"A practical extension to microfacet theory for the modeling of varying iridescence,\" ACM Transactions on Graphics (SIGGRAPH 2017 Proceedings), Volume 36, Issue 4, July 2017.";
    }

    std::string getDescription() const override
    {
        return "Reference: " + getReference();
    }

private:
    double roughness_;
    double filmThickness_;
    double filmIor_;
    double baseIor_;
    double baseExtinctionCoefficient_;
};

} // namespace lb

#endif // LIBBSDF_VARYING_IRIDESCENCE_MICROFACET_H
