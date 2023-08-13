// =================================================================== //
// Copyright (C) 2021-2023 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#include <libbsdf/Fitter/LambertianFitter.h>

#include <thread>

#include <ceres/ceres.h>

using namespace lb;

struct Cost
{
    Cost(const BrdfFitter::Sample& sample) : sample_(&sample) {}

    template <typename T>
    bool operator()(const T* const color, T* residual) const
    {
        using JetVec3 = Eigen::Matrix<T, 3, 1>;

        JetVec3 c(color[0], color[1], color[2]);
        JetVec3 value = Lambertian::compute(c);

        JetVec3 diff = BrdfFitter::toLogScale(sample_->value) - BrdfFitter::toLogScale(value);
        residual[0] = diff[0];
        residual[1] = diff[1];
        residual[2] = diff[2];

        return true;
    }

private:
    const BrdfFitter::Sample* sample_;
};

Lambertian LambertianFitter::estimateParameters(const Brdf&         brdf,
                                                int                 numSampling,
                                                const Vec3::Scalar& maxTheta)
{
    Lambertian model(Vec3(0.5, 0.5, 0.5));

    estimateParameters(&model, brdf, numSampling, maxTheta);

    return model;
}

void LambertianFitter::estimateParameters(Lambertian*         model,
                                          const Brdf&         brdf,
                                          int                 numSampling,
                                          const Vec3::Scalar& maxTheta)
{
    Data data(brdf, numSampling, maxTheta);

    ReflectanceModel::Parameters& params = model->getParameters();

    Vec3*  colorVec3 = params.at(0).getVec3();
    double color[3] = {(*colorVec3)[0], (*colorVec3)[1], (*colorVec3)[2]};

    ceres::Problem problem;

    for (auto& s : data.getSamples()) {
        Cost* cost = new Cost(s);
        ceres::CostFunction* costFunc = new ceres::AutoDiffCostFunction<Cost, 3, 3>(cost);
        problem.AddResidualBlock(costFunc, nullptr, color);
    }

    setParameterBounds(&problem, color, params.at(0));

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

    (*colorVec3)[0] = color[0];
    (*colorVec3)[1] = color[1];
    (*colorVec3)[2] = color[2];
}
