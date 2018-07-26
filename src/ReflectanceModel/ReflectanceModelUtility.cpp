// =================================================================== //
// Copyright (C) 2016-2018 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#include <libbsdf/ReflectanceModel/ReflectanceModelUtility.h>

using namespace lb;

bool reflectance_model_utility::setupTabularBrdf(const ReflectanceModel&    model,
                                                 Brdf*                      brdf,
                                                 DataType                   dataType,
                                                 float                      maxValue)
{
    using std::max;
    using std::min;

    SampleSet* ss = brdf->getSampleSet();

    ColorModel cm = ss->getColorModel();
    if (cm != RGB_MODEL &&
        cm != MONOCHROMATIC_MODEL) {
        std::cerr << "[setupTabularBrdf] Unsupported color model: " << cm << std::endl;
        return false;
    }

    for (int i0 = 0; i0 < ss->getNumAngles0(); ++i0) {
    for (int i1 = 0; i1 < ss->getNumAngles1(); ++i1) {
    for (int i2 = 0; i2 < ss->getNumAngles2(); ++i2) {
        Vec3 inDir, outDir;
        Vec3 values;
        Spectrum sp;
        #pragma omp parallel for private(inDir, outDir, values, sp)
        for (int i3 = 0; i3 < ss->getNumAngles3(); ++i3) {
            brdf->getInOutDirection(i0, i1, i2, i3, &inDir, &outDir);

            const float minZ = 0.001f;
            inDir.z() = max(inDir.z(), minZ);
            outDir.z() = max(outDir.z(), minZ);

            inDir.normalize();
            outDir.normalize();

            if (dataType == BTDF_DATA) {
                outDir.z() = -outDir.z();
            }

            values = model.getBrdfValue(inDir, outDir);
            assert(values.allFinite());

            if (cm == RGB_MODEL) {
                sp = values.asVector3f();
                sp = sp.cwiseMin(maxValue);
                ss->setSpectrum(i0, i1, i2, i3, sp);
            }
            else { // MONOCHROMATIC_MODEL
                sp.resize(1);
                sp[0] = min(values.sum() / 3.0f, maxValue);
                ss->setSpectrum(i0, i1, i2, i3, sp);
            }
        }
    }}}

    return true;
}
