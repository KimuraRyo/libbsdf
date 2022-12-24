// =================================================================== //
// Copyright (C) 2021-2022 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#include <libbsdf/Fitter/UnrealEngine4Fitter.h>

#include <thread>

#include <ceres/ceres.h>

using namespace lb;

struct Cost
{
    Cost(const BrdfFitter::Sample& sample) : sample_(&sample) {}

    template <typename T>
    bool operator()(const T* const color,
                    const T* const metallic,
                    const T* const specular,
                    const T* const roughness,
                    T*             residual) const
    {
        using JetVec3 = Eigen::Matrix<T, 3, 1>;

        const Vec3 normal(0.0, 0.0, 1.0);

        JetVec3 c(color[0], color[1], color[2]);
        JetVec3 value = UnrealEngine4::compute(sample_->inDir, sample_->outDir, normal, c,
                                               *metallic, *specular, *roughness);

        JetVec3 diff = BrdfFitter::toLogScale(sample_->value) - BrdfFitter::toLogScale(value);
        residual[0] = diff[0];
        residual[1] = diff[1];
        residual[2] = diff[2];

        return true;
    }

private:
    const BrdfFitter::Sample* sample_;
};

struct ConstMetallicCost
{
    ConstMetallicCost(const BrdfFitter::Sample& sample, const double& metallic)
        : sample_(&sample), metallic_(metallic)
    {
    }

    template <typename T>
    bool operator()(const T* const color,
                    const T* const specular,
                    const T* const roughness,
                    T*             residual) const
    {
        using JetVec3 = Eigen::Matrix<T, 3, 1>;

        const Vec3 normal(0.0, 0.0, 1.0);

        JetVec3 c(color[0], color[1], color[2]);
        JetVec3 value = UnrealEngine4::compute(sample_->inDir, sample_->outDir, normal, c,
                                               T(metallic_), *specular, *roughness);

        JetVec3 diff = BrdfFitter::toLogScale(sample_->value) - BrdfFitter::toLogScale(value);
        residual[0] = diff[0];
        residual[1] = diff[1];
        residual[2] = diff[2];

        return true;
    }

private:
    const BrdfFitter::Sample* sample_;
    double                    metallic_;
};

UnrealEngine4 UnrealEngine4Fitter::estimateParameters(const Brdf&         brdf,
                                                      int                 numSampling,
                                                      const Vec3::Scalar& maxTheta)
{
    UnrealEngine4 model(Vec3(0.5, 0.5, 0.5), 0.5f, 0.5f, 0.01f);

    estimateParameters(&model, brdf, numSampling, maxTheta);

    return model;
}

void UnrealEngine4Fitter::estimateParameters(UnrealEngine4*      model,
                                             const Brdf&         brdf,
                                             int                 numSampling,
                                             const Vec3::Scalar& maxTheta)
{
    Data data(brdf, numSampling, maxTheta);

    ReflectanceModel::Parameters& params = model->getParameters();

    UnrealEngine4 model0(*params.at(0).getVec3(),
                         0.0f,
                         *params.at(2).getFloat(),
                         *params.at(3).getFloat());

    UnrealEngine4 model1(*params.at(0).getVec3(),
                         1.0f,
                         *params.at(2).getFloat(),
                         *params.at(3).getFloat());

    Vec3*  colorVec3 = params.at(0).getVec3();
    double color[3] = {(*colorVec3)[0], (*colorVec3)[1], (*colorVec3)[2]};
    double metallic = *params.at(1).getFloat();
    double specular = *params.at(2).getFloat();
    double roughness = *params.at(3).getFloat();

    ceres::Problem problem;

    for (auto& s : data.getSamples()) {
        Cost* cost = new Cost(s);
        ceres::CostFunction* costFunc = new ceres::AutoDiffCostFunction<Cost, 3, 3, 1, 1, 1>(cost);
        problem.AddResidualBlock(costFunc, nullptr, color, &metallic, &specular, &roughness);
    }

    setParameterBounds(&problem, color, params.at(0));
    setParameterBounds(&problem, &metallic, params.at(1));
    setParameterBounds(&problem, &specular, params.at(2));
    setParameterBounds(&problem, &roughness, params.at(3));

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
    *params.at(1).getFloat() = static_cast<float>(metallic);
    *params.at(2).getFloat() = static_cast<float>(specular);
    *params.at(3).getFloat() = static_cast<float>(roughness);

    // Compare the results and select the best parameters.
    estimateWithoutMetallic(&model0, data);
    estimateWithoutMetallic(&model1, data);

    Scalar err = computeError(*model, data);
    Scalar err0 = computeError(model0, data);
    Scalar err1 = computeError(model1, data);

    if (err > err0 || err > err1) {
        ReflectanceModel::Parameters* newParams;

        if (err0 < err1) {
            newParams = &model0.getParameters();
        }
        else {
            newParams = &model1.getParameters();
        }

        *params.at(0).getVec3() = *newParams->at(0).getVec3();
        *params.at(1).getFloat() = *newParams->at(1).getFloat();
        *params.at(2).getFloat() = *newParams->at(2).getFloat();
        *params.at(3).getFloat() = *newParams->at(3).getFloat();
    }
}

void UnrealEngine4Fitter::estimateWithoutMetallic(UnrealEngine4* model, const Data& data)
{
    ReflectanceModel::Parameters& params = model->getParameters();

    Vec3*  colorVec3 = params.at(0).getVec3();
    double color[3] = {(*colorVec3)[0], (*colorVec3)[1], (*colorVec3)[2]};
    double metallic = *params.at(1).getFloat();
    double specular = *params.at(2).getFloat();
    double roughness = *params.at(3).getFloat();

    ceres::Problem problem;

    for (auto& s : data.getSamples()) {
        ConstMetallicCost* cost = new ConstMetallicCost(s, metallic);
        ceres::CostFunction* costFunc =
            new ceres::AutoDiffCostFunction<ConstMetallicCost, 3, 3, 1, 1>(cost);
        problem.AddResidualBlock(costFunc, nullptr, color, &specular, &roughness);
    }

    setParameterBounds(&problem, color, params.at(0));
    setParameterBounds(&problem, &specular, params.at(2));
    setParameterBounds(&problem, &roughness, params.at(3));

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
    *params.at(2).getFloat() = static_cast<float>(specular);
    *params.at(3).getFloat() = static_cast<float>(roughness);
}
