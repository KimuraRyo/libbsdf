// =================================================================== //
// Copyright (C) 2016 Kimura Ryo                                       //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#include <libbsdf/ReflectanceModel/ReflectanceModelUtility.h>

using namespace lb;

bool reflectance_model_utility::setupTabularBrdf(const ReflectanceModel& model,
                                                 Brdf*                   brdf)
{
    SampleSet* ss = brdf->getSampleSet();

    if (ss->getNumWavelengths() != 3) {
        std::cerr << "[setupBrdf] Currently, BRDF not using RGB mode is not supported." << std::endl;
        return false;
    }

    for (int i0 = 0; i0 < ss->getNumAngles0(); ++i0) {
    for (int i1 = 0; i1 < ss->getNumAngles1(); ++i1) {
    for (int i2 = 0; i2 < ss->getNumAngles2(); ++i2) {
    for (int i3 = 0; i3 < ss->getNumAngles3(); ++i3) {
        Vec3 inDir, outDir;
        brdf->getInOutDirection(i0, i1, i2, i3, &inDir, &outDir);

        const float minZ = 0.001f;
        inDir.z() = std::max(inDir.z(), minZ);
        outDir.z() = std::max(outDir.z(), minZ);

        inDir.normalize();
        outDir.normalize();

        Vec3 values = model.getBrdfValue(inDir, outDir);
        assert(values.allFinite());

        Spectrum sp = values.asVector3f();
        const float maxBrdfVal = 10000.0f;
        sp = sp.cwiseMin(maxBrdfVal);
        ss->setSpectrum(i0, i1, i2, i3, sp);
    }}}}

    return true;
}
