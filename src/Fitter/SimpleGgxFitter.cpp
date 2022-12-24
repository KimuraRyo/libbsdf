// =================================================================== //
// Copyright (C) 2022 Kimura Ryo                                       //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#include <libbsdf/Fitter/SimpleGgxFitter.h>

#include <thread>

#include <ceres/ceres.h>

using namespace lb;

struct Cost
{
    Cost(const BrdfFitter::Sample& sample) : sample_(&sample) {}

    template <typename T>
    bool operator()(const T* const color, const T* const roughness, T* residual) const
    {
        using JetVec3 = Eigen::Matrix<T, 3, 1>;

        JetVec3 inDir(T(sample_->inDir[0]), T(sample_->inDir[1]), T(sample_->inDir[2]));
        JetVec3 outDir(T(sample_->outDir[0]), T(sample_->outDir[1]), T(sample_->outDir[2]));
        const JetVec3 normal(T(0), T(0), T(1));

        JetVec3 c(color[0], color[1], color[2]);
        JetVec3 value = SimpleGgx::compute(inDir, outDir, normal, c, *roughness);

        JetVec3 diff = BrdfFitter::toLogScale(sample_->value) - BrdfFitter::toLogScale(value);
        residual[0] = diff[0];
        residual[1] = diff[1];
        residual[2] = diff[2];

        return true;
    }

private:
    const BrdfFitter::Sample* sample_;
};

SimpleGgx
SimpleGgxFitter::estimateParameters(const Brdf& brdf, int numSampling, const Vec3::Scalar& maxTheta)
{
    SimpleGgx model(Vec3(0.5, 0.5, 0.5), 0.01f);

    estimateParameters(&model, brdf, numSampling, maxTheta);

    return model;
}

void SimpleGgxFitter::estimateParameters(SimpleGgx*          model,
                                         const Brdf&         brdf,
                                         int                 numSampling,
                                         const Vec3::Scalar& maxTheta)
{
    Data data(brdf, numSampling, maxTheta);

    ReflectanceModel::Parameters& params = model->getParameters();

    Vec3*  colorVec3 = params.at(0).getVec3();
    double color[3] = {(*colorVec3)[0], (*colorVec3)[1], (*colorVec3)[2]};
    double roughness = *params.at(1).getFloat();

    ceres::Problem problem;

    for (auto& s : data.getSamples()) {
        Cost* cost = new Cost(s);
        ceres::CostFunction* costFunc = new ceres::AutoDiffCostFunction<Cost, 3, 3, 1>(cost);
        problem.AddResidualBlock(costFunc, nullptr, color, &roughness);
    }

    setParameterBounds(&problem, color, params.at(0));
    setParameterBounds(&problem, &roughness, params.at(1));

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
    *params.at(1).getFloat() = static_cast<float>(roughness);
}
