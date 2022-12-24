// =================================================================== //
// Copyright (C) 2021-2022 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#include <libbsdf/Fitter/BrdfFitter.h>

#include <ceres/ceres.h>

#include <libbsdf/Brdf/SpecularCoordinatesBrdf.h>
#include <libbsdf/Brdf/SphericalCoordinatesBrdf.h>
#include <libbsdf/Common/SpectrumUtility.h>
#include <libbsdf/Common/Utility.h>
#include <libbsdf/Common/Xorshift.h>

using namespace lb;

BrdfFitter::Data::Data(const Brdf& brdf, int numSampling, const Vec3::Scalar& maxTheta)
{
    if (maxTheta < toRadian(Vec3::Scalar(10))) {
        lbWarn << "[BrdfFitter::Data::Data] maxTheta is too small. maxTheta: " << maxTheta;
    }

    const SampleSet* ss = brdf.getSampleSet();

    const Vec3::Scalar threshold = std::max(std::cos(maxTheta), Vec3::Scalar(0));

    if (numSampling <= 0 || maxTheta < toRadian(Vec3::Scalar(EPSILON_F))) {
        for (int i0 = 0; i0 < ss->getNumAngles0(); ++i0) {
        for (int i1 = 0; i1 < ss->getNumAngles1(); ++i1) {
        for (int i2 = 0; i2 < ss->getNumAngles2(); ++i2) {
        for (int i3 = 0; i3 < ss->getNumAngles3(); ++i3) {
            Vec3 inDir, outDir;
            brdf.toXyz(ss->getAngle0(i0),
                       ss->getAngle1(i1),
                       ss->getAngle2(i2),
                       ss->getAngle3(i3),
                       &inDir, &outDir);

            if (inDir[2] < threshold || outDir[2] < threshold) continue;

            Sample sample;
            sample.inDir = inDir;
            sample.outDir = outDir;
            sample.value = SpectrumUtility::toSrgb<Vec3, SampleSet>(
                brdf.getSpectrum(sample.inDir, sample.outDir), *ss);

            samples_.push_back(sample);
        }}}}
    }
    else {
        int count = 0;
        Xorshift rnd;

        bool inDirIsParam = dynamic_cast<const SpecularCoordinatesBrdf*>(&brdf) ||
                            dynamic_cast<const SphericalCoordinatesBrdf*>(&brdf);

        while (count < numSampling) {
            Sample sample;

            // If the type of BRDF has incident directions as parameters, they are taken as sample points for fitting.
            if (inDirIsParam) {
                int inThIndex = Xorshift::random() % ss->getNumAngles0();

                int inPhIndex;
                if (ss->isIsotropic()) {
                    inPhIndex = 0;
                }
                else {
                    inPhIndex = Xorshift::random() % (ss->getNumAngles1() - 1);
                }

                sample.inDir = SphericalCoordinateSystem::toXyz(ss->getAngle0(inThIndex),
                                                                ss->getAngle1(inPhIndex));
            }
            else {
                sample.inDir = rnd.nextOnHemisphere<Vec3>();
            }

            sample.outDir = rnd.nextOnHemisphere<Vec3>();

            if (sample.inDir[2] < threshold || sample.outDir[2] < threshold) continue;

            sample.value = SpectrumUtility::toSrgb<Vec3, SampleSet>(
                brdf.getSpectrum(sample.inDir, sample.outDir), *ss);

            samples_.push_back(sample);
            ++count;
        }
    }
}

bool BrdfFitter::setParameterBounds(ceres::Problem*                    problem,
                                    double*                            parameter,
                                    const ReflectanceModel::Parameter& bounds)
{
    switch (bounds.getType()) {
        case ReflectanceModel::Parameter::FLOAT_PARAMETER:
            problem->SetParameterLowerBound(parameter, 0, *bounds.getMinFloat());
            problem->SetParameterUpperBound(parameter, 0, *bounds.getMaxFloat());
            break;
        case ReflectanceModel::Parameter::VEC3_PARAMETER: {
            const Vec3* minVec3 = bounds.getMinVec3();
            const Vec3* maxVec3 = bounds.getMaxVec3();
            problem->SetParameterLowerBound(parameter, 0, (*minVec3)[0]);
            problem->SetParameterLowerBound(parameter, 1, (*minVec3)[1]);
            problem->SetParameterLowerBound(parameter, 2, (*minVec3)[2]);
            problem->SetParameterUpperBound(parameter, 0, (*maxVec3)[0]);
            problem->SetParameterUpperBound(parameter, 1, (*maxVec3)[1]);
            problem->SetParameterUpperBound(parameter, 2, (*maxVec3)[2]);
            break;
        }
        case ReflectanceModel::Parameter::INT_PARAMETER:
            break;
        default:
            lbError << "[BrdfFitter::setParameterBounds] Invalid parameter type: " << bounds.getType();
            return false;
    }

    return true;
}

//Vec3::Scalar BrdfFitter::computeRmse(const ReflectanceModel& model, const Data& data)
//{
//    Vec3d sum = Vec3d::Zero();
//    for (auto& s : data.getSamples()) {
//        Vec3 diff = model.getBrdfValue(s.inDir, s.outDir) - s.value;
//        sum += diff.cwiseProduct(diff);
//    }
//
//    Vec3d avg = sum / data.getSamples().size();
//    Vec3d rmse = avg.cwiseSqrt();
//    return rmse.sum();
//}

Vec3::Scalar BrdfFitter::computeError(const ReflectanceModel& model, const Data& data)
{
    Vec3d sum = Vec3d::Zero();
    for (auto& s : data.getSamples()) {
        Vec3 diff = toLogScale(model.getBrdfValue(s.inDir, s.outDir)) - toLogScale(s.value);
        sum += diff.cwiseAbs();
    }

    Vec3d avg = sum / data.getSamples().size();

    return avg.sum();
}
