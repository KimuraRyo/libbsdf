// =================================================================== //
// Copyright (C) 2014-2015 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#include <libbsdf/Brdf/Processor.h>
#include <libbsdf/Common/SphericalCoordinateSystem.h>

using namespace lb;

void lb::divideByCosineOutTheta(Brdf* brdf)
{
    SampleSet* ss = brdf->getSampleSet();

    for (int i0 = 0; i0 < ss->getNumAngles0(); ++i0) {
    for (int i1 = 0; i1 < ss->getNumAngles1(); ++i1) {
    for (int i2 = 0; i2 < ss->getNumAngles2(); ++i2) {
    for (int i3 = 0; i3 < ss->getNumAngles3(); ++i3) {
        Vec3 inDir, outDir;
        brdf->getInOutDirection(i0, i1, i2, i3, &inDir, &outDir);
        float cosOutTheta = outDir.dot(Vec3(0.0, 0.0, 1.0));

        lb::Spectrum& sp = ss->getSpectrum(i0, i1, i2, i3);

        // Copy the spectrum if the Z-component of the outgoing direction is zero or negative.
        bool isZeroOrLess = (cosOutTheta <= std::cos(SphericalCoordinateSystem::MAX_ANGLE2));
        if (isZeroOrLess && i2 > 0) {
            // Assume i2 is the index of the polar angle related to outgoing directions.
            brdf->getInOutDirection(i0, i1, i2 - 1, i3, &inDir, &outDir);
            sp = ss->getSpectrum(i0, i1, i2 - 1, i3);
            cosOutTheta = outDir.dot(Vec3(0.0, 0.0, 1.0));
        }

        sp /= cosOutTheta;
    }}}}
}
