// =================================================================== //
// Copyright (C) 2016-2023 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#include <libbsdf/Brdf/SpecularCoordinatesRandomSampleSet.h>

using namespace lb;

void SpecularCoordinatesRandomSampleSet::setupBrdf(SpecularCoordinatesBrdf* brdf,
                                                   double                   weight0,
                                                   double                   weight1,
                                                   double                   weight2,
                                                   double                   weight3)
{
    for (int inThIndex = 0; inThIndex < brdf->getNumInTheta();   ++inThIndex) {
    for (int inPhIndex = 0; inPhIndex < brdf->getNumInPhi();     ++inPhIndex) {
    for (int spThIndex = 0; spThIndex < brdf->getNumSpecTheta(); ++spThIndex) {
        double coeff3;
        if (brdf->getSampleSet()->isIsotropic()) {
            // Compute a weight coefficient for specular azimuthal angles along specular polar angle.
            coeff3 = brdf->getSpecTheta(spThIndex) / SpecularCoordinateSystem::MAX_ANGLE2;
        }
        else {
            coeff3 = 1;
        }

        RandomSampleSet::AngleList angles;
        SampleMap::iterator it;
        double w3;
        Spectrum sp;
        #pragma omp parallel for private(angles, it, w3, sp)
        for (int spPhIndex = 0; spPhIndex < brdf->getNumSpecPhi(); ++spPhIndex) {
            angles.resize(4);
            angles.at(0) = brdf->getInTheta(inThIndex);
            angles.at(1) = brdf->getInPhi(inPhIndex);
            angles.at(2) = brdf->getSpecTheta(spThIndex);
            angles.at(3) = brdf->getSpecPhi(spPhIndex);

            it = sampleMap_.find(angles);
            if (it != sampleMap_.end()) {
                brdf->setSpectrum(inThIndex, inPhIndex, spThIndex, spPhIndex, it->second);
            }
            else {
                // Modify a weight coefficient for specular azimuthal angles.
                w3 = hermiteInterpolation3(0.0, weight3, coeff3);
                sp = estimateSpectrum<SpecularCoordinateSystem>(angles, weight0, weight1, weight2, w3);
                brdf->setSpectrum(inThIndex, inPhIndex, spThIndex, spPhIndex, sp);
            }
        }
    }}}

    brdf->getSampleSet()->updateAngleAttributes();
}
