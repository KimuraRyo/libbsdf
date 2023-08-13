// =================================================================== //
// Copyright (C) 2018-2023 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#ifndef LIBBSDF_MULTIPLE_SCATTERING_SMITH_H
#define LIBBSDF_MULTIPLE_SCATTERING_SMITH_H

#include <libbsdf/Common/Global.h>
#include <libbsdf/ReflectanceModel/ReflectanceModel.h>

namespace lb {

/*! Multiple scattering Smith reflectance model. */
class MultipleScatteringSmith : public ReflectanceModel
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    /*! Types of material. */
    enum MaterialType
    {
        CONDUCTOR_MATERIAL = 0,
        DIELECTRIC_MATERIAL = 1,
        DIFFUSE_MATERIAL = 2
    };

    /*! Types of height distribution. */
    enum HeightType
    {
        UNIFORM_HEIGHT = 0,
        GAUSSIAN_HEIGHT = 1
    };

    /*! Types of slope distribution. */
    enum SlopeType
    {
        BECKMANN_SLOPE = 0,
        GGX_SLOPE = 1
    };

    MultipleScatteringSmith(const Vec3& color,
                            double      roughnessX,
                            double      roughnessY,
                            double      refractiveIndex,
                            int         materialType,
                            int         heightType,
                            int         slopeType,
                            int         numIterations)
        : color_(color),
          roughnessX_(roughnessX),
          roughnessY_(roughnessY),
          refractiveIndex_(refractiveIndex),
          materialType_(static_cast<MaterialType>(materialType)),
          heightType_(static_cast<HeightType>(heightType)),
          slopeType_(static_cast<SlopeType>(slopeType)),
          numIterations_(numIterations)
    {
        parameters_.push_back(Parameter("Color", &color_));
        parameters_.push_back(Parameter("Roughness X", &roughnessX_, 0.01, 100.0));
        parameters_.push_back(Parameter("Roughness Y", &roughnessY_, 0.01, 100.0));
        parameters_.push_back(Parameter("Refractive index", &refractiveIndex_, 0.01, 100.0,
                                        "Only valid for dielectric material"));
        parameters_.push_back(Parameter("Material type", &materialType_, 0, 2,
                                        "0: Conductor\n1: Dielectric\n2: Diffuse"));
        parameters_.push_back(
            Parameter("Height type", &heightType_, 0, 1, "0: Uniform\n1: Gaussian"));
        parameters_.push_back(Parameter("Slope type", &slopeType_, 0, 1, "0: Beckmann\n1: GGX"));
        parameters_.push_back(Parameter("Number of iterations", &numIterations_));
    }

    /*! Evaluates a BSDF value with iterations. */
    static Vec3 compute(const Vec3&  L,
                        const Vec3&  V,
                        const Vec3&  color,
                        double       roughnessX,
                        double       roughnessY,
                        double       refractiveIndex,
                        MaterialType materialType,
                        HeightType   heightType,
                        SlopeType    slopeType,
                        int          numIterations);

    Vec3 getValue(const Vec3& inDir, const Vec3& outDir) const
    {
        return compute(inDir, outDir, color_, roughnessX_, roughnessY_, refractiveIndex_,
                       static_cast<MaterialType>(materialType_),
                       static_cast<HeightType>(heightType_), static_cast<SlopeType>(slopeType_),
                       numIterations_);
    }

    bool isIsotropic() const { return false; }

    std::string getName() const { return "Multiple scattering Smith"; }

    static std::string getReference()
    {
        return "Eric Heitz, Johannes Hanika, Eugene d'Eon, and Carsten Dachsbacher, \"Multiple-scattering microfacet BSDFs with the Smith model,\" ACM Transactions on Graphics (SIGGRAPH 2016 Proceedings), Volume 35, Issue 4, July 2016.";
    }

    std::string getDescription() const
    {
        return "Reference: " + getReference();
    }

private:
    Vec3   color_;
    double roughnessX_;
    double roughnessY_;
    double refractiveIndex_;
    int    numIterations_;
    int    materialType_;
    int    heightType_;
    int    slopeType_;
};

} // namespace lb

#endif // LIBBSDF_MULTIPLE_SCATTERING_SMITH_H
