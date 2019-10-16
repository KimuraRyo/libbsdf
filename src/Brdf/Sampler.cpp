// =================================================================== //
// Copyright (C) 2016-2019 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#include <libbsdf/Brdf/Sampler.h>

using namespace lb;

bool Sampler::isIsotropic(const SampleSet& samples)
{
    return samples.isIsotropic();
}

bool Sampler::isIsotropic(const SampleSet2D& ss2)
{
    return ss2.isIsotropic();
}

const SampleSet* Sampler::getSampleSet(const Brdf& brdf)
{
    return brdf.getSampleSet();
}

void Sampler::fromXyz(const Brdf& brdf,
                      const Vec3& inDir, const Vec3& outDir,
                      float* angle0, float* angle2, float* angle3)
{
    brdf.fromXyz(inDir, outDir, angle0, angle2, angle3);
}

void Sampler::fromXyz(const Brdf& brdf,
                      const Vec3& inDir, const Vec3& outDir,
                      float* angle0, float* angle1, float* angle2, float* angle3)
{
    brdf.fromXyz(inDir, outDir, angle0, angle1, angle2, angle3);
}
