// =================================================================== //
// Copyright (C) 2016-2019 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#include <libbsdf/ReflectanceModel/ReflectanceModelUtility.h>

#include <cassert>

using namespace lb;

bool reflectance_model_utility::setupTabularBrdf(const ReflectanceModel&    model,
                                                 Brdf*                      brdf,
                                                 DataType                   dataType,
                                                 float                      maxValue)
{
    using std::abs;
    using std::max;
    using std::min;

    SampleSet* ss = brdf->getSampleSet();

    ColorModel cm = ss->getColorModel();
    if (cm != RGB_MODEL &&
        cm != MONOCHROMATIC_MODEL) {
        lbError << "[reflectance_model_utility::setupTabularBrdf] Unsupported color model: " << cm;
        return false;
    }

    Vec3 inDir, outDir;
    Vec3 values;
    Spectrum sp;
    int i0, i1, i3;
    #pragma omp parallel for private(inDir, outDir, values, sp, i0, i1, i3) schedule(dynamic)
    for (int i2 = 0; i2 < ss->getNumAngles2(); ++i2) {
        for (i0 = 0; i0 < ss->getNumAngles0(); ++i0) {
        for (i1 = 0; i1 < ss->getNumAngles1(); ++i1) {
        for (i3 = 0; i3 < ss->getNumAngles3(); ++i3) {
            brdf->getInOutDirection(i0, i1, i2, i3, &inDir, &outDir);

            // Adjust horizontal and downward directions.
            const Vec3::Scalar epsilon = Vec3::Scalar(0.001);
            inDir.z()  = max(inDir.z(),  epsilon);
            outDir.z() = max(outDir.z(), epsilon);

            // Adjust a downward outgoing direction along the Z-axis.
            if (abs(outDir.x()) <= epsilon &&
                abs(outDir.y()) <= epsilon &&
                outDir.z() <= epsilon) {
                outDir.x() = 1;
            }

            inDir.normalize();
            outDir.normalize();

            if (dataType == BTDF_DATA) {
                outDir.z() = -outDir.z();
            }

            values = model.getBrdfValue(inDir, outDir);
            assert(values.allFinite());

            if (cm == RGB_MODEL) {
                sp = values.cast<Spectrum::Scalar>();
                sp = sp.cwiseMin(maxValue);
            }
            else { // MONOCHROMATIC_MODEL
                sp.resize(1);
                sp[0] = min(static_cast<float>(values.sum()) / 3.0f, maxValue);
            }
            ss->setSpectrum(i0, i1, i2, i3, sp);
        }}}
    }

    return true;
}
