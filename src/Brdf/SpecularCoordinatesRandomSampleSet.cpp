// =================================================================== //
// Copyright (C) 2016 Kimura Ryo                                       //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#include <libbsdf/Brdf/SpecularCoordinatesRandomSampleSet.h>

using namespace lb;

void SpecularCoordinatesRandomSampleSet::setupBrdf(SpecularCoordinatesBrdf* brdf,
                                                   float                    weight0,
                                                   float                    weight1,
                                                   float                    weight2,
                                                   float                    weight3)
{
    for (int inThIndex = 0; inThIndex < brdf->getNumInTheta();   ++inThIndex) {
    for (int inPhIndex = 0; inPhIndex < brdf->getNumInPhi();     ++inPhIndex) {
    for (int spThIndex = 0; spThIndex < brdf->getNumSpecTheta(); ++spThIndex) {
        RandomSampleSet::AngleList angles;
        SampleMap::iterator it;
        float w3;
        #pragma omp parallel for private(angles, it, w3)
        for (int spPhIndex = 0; spPhIndex < brdf->getNumSpecPhi(); ++spPhIndex) {
            angles.resize(4);
            angles.at(0) = brdf->getInTheta(inThIndex);
            angles.at(1) = brdf->getInPhi(inPhIndex);
            angles.at(2) = brdf->getSpecTheta(spThIndex);
            angles.at(3) = brdf->getSpecPhi(spPhIndex);

            it = sampleMap_.find(angles);
            if (it != sampleMap_.end()) {
                brdf->setSpectrum(inThIndex, inPhIndex, spThIndex, spPhIndex,
                                  it->second);
            }
            else {
                // Modify a weight coefficient.
                w3 = weight3 * hermiteInterpolation3(1.0f, 1.0f / weight3, brdf->getSpecTheta(spThIndex) / PI_2_F);
                Spectrum sp = estimateSpectrum<SpecularCoordinateSystem>(angles, weight0, weight1, weight2, w3);
                brdf->setSpectrum(inThIndex, inPhIndex, spThIndex, spPhIndex, sp);
            }
        }
    }}}

    brdf->getSampleSet()->updateAngleAttributes();
}
