// =================================================================== //
// Copyright (C) 2014-2019 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#include <libbsdf/Brdf/Integrator.h>

#include <libbsdf/Common/Xorshift.h>

using namespace lb;

Integrator::Integrator(int numSampling) : numSampling_(numSampling)
{
    initializeOutDirs();
}

Spectrum Integrator::computeReflectance(const Brdf& brdf, const Vec3& inDir)
{
    Arrayd sumSpectrum;
    sumSpectrum.resize(brdf.getSampleSet()->getNumWavelengths());
    sumSpectrum.setZero();

    Vec3 outDir;
    Spectrum sp;
    #pragma omp parallel for private(outDir, sp)
    for (int i = 0; i < numSampling_; ++i) {
        outDir = outDirs_.col(i).cast<Vec3::Scalar>();
        sp = brdf.getSpectrum(inDir, outDir);
        sp *= static_cast<Spectrum::Scalar>(outDir.z());

        #pragma omp critical
        sumSpectrum += sp.cast<Arrayd::Scalar>();
    }

    sumSpectrum *= TAU_D / numSampling_;
    return sumSpectrum.cast<Spectrum::Scalar>();
}

Spectrum Integrator::computeReflectance(const Brdf& brdf, const Vec3& inDir, int numSampling)
{
    Arrayd sumSpectrum;
    sumSpectrum.resize(brdf.getSampleSet()->getNumWavelengths());
    sumSpectrum.setZero();

    Vec3 outDir;
    Spectrum sp;
    #pragma omp parallel for private(outDir, sp)
    for (int i = 0; i < numSampling; ++i) {
        outDir = Xorshift::randomOnHemisphere<Vec3>();
        sp = brdf.getSpectrum(inDir, outDir);
        sp *= static_cast<Spectrum::Scalar>(outDir.z());

        #pragma omp critical
        sumSpectrum += sp.cast<Arrayd::Scalar>();
    }

    sumSpectrum *= TAU_D / numSampling;
    return sumSpectrum.cast<Spectrum::Scalar>();
}

void Integrator::initializeOutDirs()
{
    outDirs_.resize(Eigen::NoChange, numSampling_);
    for (int i = 0; i < numSampling_; ++i) {
        outDirs_.col(i) = Xorshift::randomOnHemisphere<Vec3f>();
    }
}
