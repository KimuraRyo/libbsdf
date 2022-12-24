// =================================================================== //
// Copyright (C) 2021-2022 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#include <libbsdf/Fitter/GgxAnisotropicFitter.h>

#include <thread>

#include <ceres/ceres.h>

using namespace lb;

struct Cost
{
    Cost(const BrdfFitter::Sample& sample) : sample_(&sample) {}

    template <typename T>
    bool operator()(const T* const color,
                    const T* const roughnessX,
                    const T* const roughnessY,
                    const T* const refractiveIndex,
                    const T* const extinctionCoefficient,
                    T*             residual) const
    {
        using JetVec3 = Eigen::Matrix<T, 3, 1>;

        JetVec3 inDir(T(sample_->inDir[0]), T(sample_->inDir[1]), T(sample_->inDir[2]));
        JetVec3 outDir(T(sample_->outDir[0]), T(sample_->outDir[1]), T(sample_->outDir[2]));
        const JetVec3 normal(T(0), T(0), T(1));
        const JetVec3 tangent(T(1), T(0), T(0));
        const JetVec3 binormal(T(0), T(-1), T(0));

        JetVec3 c(color[0], color[1], color[2]);
        JetVec3 value =
            GgxAnisotropic::compute(inDir, outDir, normal, tangent, binormal, c, *roughnessX,
                                    *roughnessY, *refractiveIndex, *extinctionCoefficient);

        JetVec3 diff = BrdfFitter::toLogScale(sample_->value) - BrdfFitter::toLogScale(value);
        residual[0] = diff[0];
        residual[1] = diff[1];
        residual[2] = diff[2];

        return true;
    }

private:
    const BrdfFitter::Sample* sample_;
};

GgxAnisotropic GgxAnisotropicFitter::estimateParameters(const Brdf&         brdf,
                                                        int                 numSampling,
                                                        const Vec3::Scalar& maxTheta)
{
    GgxAnisotropic model(Vec3(0.5, 0.5, 0.5), 0.01f, 0.01f, 1.5f, 0.0f);

    estimateParameters(&model, brdf, numSampling, maxTheta);

    return model;
}

void GgxAnisotropicFitter::estimateParameters(GgxAnisotropic*     model,
                                              const Brdf&         brdf,
                                              int                 numSampling,
                                              const Vec3::Scalar& maxTheta)
{
    Data data(brdf, numSampling, maxTheta);

    ReflectanceModel::Parameters& params = model->getParameters();

    Vec3*  colorVec3 = params.at(0).getVec3();
    double color[3] = {(*colorVec3)[0], (*colorVec3)[1], (*colorVec3)[2]};
    double roughnessX = *params.at(1).getFloat();
    double roughnessY = *params.at(2).getFloat();
    double refractiveIndex = *params.at(3).getFloat();
    double extinctionCoefficient = *params.at(4).getFloat();

    ceres::Problem problem;

    for (auto& s : data.getSamples()) {
        Cost* cost = new Cost(s);
        ceres::CostFunction* costFunc =
            new ceres::AutoDiffCostFunction<Cost, 3, 3, 1, 1, 1, 1>(cost);
        problem.AddResidualBlock(costFunc, nullptr, color, &roughnessX, &roughnessY,
                                 &refractiveIndex, &extinctionCoefficient);
    }

    setParameterBounds(&problem, color, params.at(0));
    setParameterBounds(&problem, &roughnessX, params.at(1));
    setParameterBounds(&problem, &roughnessY, params.at(2));
    setParameterBounds(&problem, &refractiveIndex, params.at(3));
    setParameterBounds(&problem, &extinctionCoefficient, params.at(4));

    ceres::Solver::Options options;
    options.linear_solver_type = ceres::DENSE_QR;
    options.minimizer_progress_to_stdout = true;
    options.num_threads = std::thread::hardware_concurrency();

    if (Log::getNotificationLevel() > Log::Level::INFO_MSG) {
        options.logging_type = ceres::SILENT;
    }

    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);
    lbInfo << summary.FullReport();

    using Scalar = typename Vec3::Scalar;

    (*colorVec3)[0] = static_cast<Scalar>(color[0]);
    (*colorVec3)[1] = static_cast<Scalar>(color[1]);
    (*colorVec3)[2] = static_cast<Scalar>(color[2]);
    *params.at(1).getFloat() = static_cast<float>(roughnessX);
    *params.at(2).getFloat() = static_cast<float>(roughnessY);
    *params.at(3).getFloat() = static_cast<float>(refractiveIndex);
    *params.at(4).getFloat() = static_cast<float>(extinctionCoefficient);
}
