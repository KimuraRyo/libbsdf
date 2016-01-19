// =================================================================== //
// Copyright (C) 2015 Kimura Ryo                                       //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#include <libbsdf/Brdf/RandomSampleSet.h>

#include <libbsdf/Common/SpecularCoordinateSystem.h>

using namespace lb;

const Spectrum& RandomSampleSet::getSpectrumOfNearestSample(const AngleList& angles,
                                                            bool             useReciprocity) const
{
    using std::acos;

    const Spectrum* sp;

    Vec3 inDir, outDir;
    SphericalCoordinateSystem::toXyz(angles.at(0), angles.at(1),
                                     angles.at(2), angles.at(3),
                                     &inDir, &outDir);

    float minAngleDiff = std::numeric_limits<float>::max();

    for (auto it = sampleMap_.begin(); it != sampleMap_.end(); ++it) {
        const AngleList& sampleAngles = it->first;
        Vec3 sampleInDir, sampleOutDir;
        SphericalCoordinateSystem::toXyz(sampleAngles.at(0), sampleAngles.at(1),
                                         sampleAngles.at(2), sampleAngles.at(3),
                                         &sampleInDir, &sampleOutDir);

        float inAngle  = acos(inDir.dot(sampleInDir));
        float outAngle = acos(outDir.dot(sampleOutDir));
        float angleDiff = inAngle + outAngle;

        if (useReciprocity) {
            float inOutAngle = acos(inDir.dot(sampleOutDir));
            float outInAngle = acos(outDir.dot(sampleInDir));
            float angleUsingReciprocity = inOutAngle + outInAngle;
            if (angleUsingReciprocity < angleDiff) {
                angleDiff = angleUsingReciprocity;
            }
        }

        if (angleDiff == 0.0f) {
            return it->second;
        }
        else if (angleDiff < minAngleDiff) {
            minAngleDiff = angleDiff;
            sp = &it->second;
        }
    }

    return *sp;
}

void RandomSampleSet::setupBrdf(SphericalCoordinatesBrdf* brdf)
{
    for (int inThIndex  = 0; inThIndex  < brdf->getNumInTheta();  ++inThIndex)  {
    for (int inPhIndex  = 0; inPhIndex  < brdf->getNumInPhi();    ++inPhIndex)  {
    for (int outThIndex = 0; outThIndex < brdf->getNumOutTheta(); ++outThIndex) {
        RandomSampleSet::AngleList angles;
        SampleMap::iterator it;
        #pragma omp parallel for private(angles, it)
        for (int outPhIndex = 0; outPhIndex < brdf->getNumOutPhi(); ++outPhIndex) {
            angles.resize(4);
            angles.at(0) = brdf->getInTheta(inThIndex);
            angles.at(1) = brdf->getInPhi(inPhIndex);
            angles.at(2) = brdf->getOutTheta(outThIndex);
            angles.at(3) = brdf->getOutPhi(outPhIndex);

            it = sampleMap_.find(angles);
            if (it != sampleMap_.end()) {
                brdf->setSpectrum(inThIndex, inPhIndex, outThIndex, outPhIndex,
                                  it->second);
            }
            else {
                brdf->setSpectrum(inThIndex, inPhIndex, outThIndex, outPhIndex,
                                  getSpectrumOfNearestSample(angles, false));
                //brdf->setSpectrum(inThIndex, inPhIndex, outThIndex, outPhIndex,
                //                  estimateSpectrum<SpecularCoordinateSystem>(angles));
            }
        }
    }}}

    brdf->getSampleSet()->updateAngleAttributes();
}
