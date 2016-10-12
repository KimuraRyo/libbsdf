// =================================================================== //
// Copyright (C) 2016 Kimura Ryo                                       //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#include <libbsdf/Brdf/SphericalCoordinatesRandomSampleSet.h>

using namespace lb;

void SphericalCoordinatesRandomSampleSet::setupBrdf(SphericalCoordinatesBrdf* brdf)
{
    for (int inThIndex = 0;  inThIndex < brdf->getNumInTheta();   ++inThIndex)  {
    for (int inPhIndex = 0;  inPhIndex < brdf->getNumInPhi();     ++inPhIndex)  {
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
                                  findSpectrumOfNearestSample(angles, false));
                //brdf->setSpectrum(inThIndex, inPhIndex, outThIndex, outPhIndex,
                //                  estimateSpectrum<CoordSysT>(angles));
            }
        }
    }}}

    brdf->getSampleSet()->updateAngleAttributes();
}
