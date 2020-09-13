// =================================================================== //
// Copyright (C) 2016-2020 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#include <libbsdf/ReflectanceModel/ReflectanceModelUtility.h>

#include <cassert>
#include <set>

#include <libbsdf/Brdf/SpecularCoordinatesBrdf.h>

using namespace lb;

bool ReflectanceModelUtility::setupBrdf(const ReflectanceModel& model,
                                        Brdf*                   brdf,
                                        DataType                dataType)
{
    SampleSet* ss = brdf->getSampleSet();

    ColorModel cm = ss->getColorModel();
    if (cm != RGB_MODEL &&
        cm != MONOCHROMATIC_MODEL) {
        lbError << "[ReflectanceModelUtility::setupBrdf] Unsupported color model: " << cm;
        return false;
    }

    Spectrum sp;
    int i0, i1, i3;
    #pragma omp parallel for private(sp, i0, i1, i3) schedule(dynamic)
    for (int i2 = 0; i2 < ss->getNumAngles2(); ++i2) {
        for (i0 = 0; i0 < ss->getNumAngles0(); ++i0) {
        for (i1 = 0; i1 < ss->getNumAngles1(); ++i1) {
        for (i3 = 0; i3 < ss->getNumAngles3(); ++i3) {
            sp = getSpectrum(model, *brdf, dataType,
                             ss->getAngle0(i0),
                             ss->getAngle1(i1),
                             ss->getAngle2(i2),
                             ss->getAngle3(i3));

            ss->setSpectrum(i0, i1, i2, i3, sp);
        }}}
    }

    return true;
}

bool ReflectanceModelUtility::setupBrdf(const ReflectanceModel& model,
                                        Brdf*                   brdf,
                                        int                     numAngles0,
                                        int                     numAngles1,
                                        int                     numAngles2,
                                        int                     numAngles3,
                                        DataType                dataType,
                                        float                   ior)
{
    SampleSet* ss = brdf->getSampleSet();

    ColorModel cm = ss->getColorModel();
    if (cm != RGB_MODEL &&
        cm != MONOCHROMATIC_MODEL) {
        lbError << "[ReflectanceModelUtility::setupBrdf] Unsupported color model: " << cm;
        return false;
    }

    Arrayf& angles0 = ss->getAngles0();
    Arrayf& angles1 = ss->getAngles1();
    Arrayf& angles2 = ss->getAngles2();
    Arrayf& angles3 = ss->getAngles3();

    bool filled0 = (angles0.size() >= numAngles0);
    bool filled1 = (angles1.size() >= numAngles1);
    bool filled2 = (angles2.size() >= numAngles2);
    bool filled3 = (angles3.size() >= numAngles3);

    while (!(filled0 && filled1 && filled2 && filled3)) {
        if (!filled0) {
            filled0 = insertAngle0(model, brdf, numAngles0, dataType);

            auto specBrdf = dynamic_cast<SpecularCoordinatesBrdf*>(brdf);
            if (specBrdf &&
                specBrdf->getNumSpecularOffsets() != 0) {
                specBrdf->setupSpecularOffsets(ior);
            }
        }

        if (!filled1) {
            filled1 = insertAngle1(model, brdf, numAngles1, dataType);
        }

        if (!filled2) {
            filled2 = insertAngle2(model, brdf, numAngles2, dataType);
        }

        if (!filled3) {
            filled3 = insertAngle3(model, brdf, numAngles3, dataType);
        }
    }

    ss->resizeAngles(static_cast<int>(angles0.size()),
                     static_cast<int>(angles1.size()),
                     static_cast<int>(angles2.size()),
                     static_cast<int>(angles3.size()));
    ss->updateAngleAttributes();

    setupBrdf(model, brdf, dataType);

    return true;
}

Spectrum ReflectanceModelUtility::getSpectrum(const ReflectanceModel&   model,
                                              const Brdf&               brdf,
                                              DataType                  dataType,
                                              float                     angle0,
                                              float                     angle1,
                                              float                     angle2,
                                              float                     angle3)
{
    using std::abs;
    using std::max;

    Vec3 inDir, outDir;
    brdf.toXyz(angle0, angle1, angle2, angle3, &inDir, &outDir);

    // Adjust horizontal and downward directions.
    const Vec3::Scalar epsilon = Vec3::Scalar(EPSILON_F);
    inDir.z()  = max(inDir.z(), Vec3::Scalar(0.001));
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

    Vec3 values = model.getBrdfValue(inDir, outDir);
    assert(values.allFinite());

    ColorModel cm = brdf.getSampleSet()->getColorModel();
    if (cm == RGB_MODEL) {
        return values.cast<Spectrum::Scalar>();
    }
    else { // MONOCHROMATIC_MODEL
        Spectrum sp(1);
        sp[0] = static_cast<float>(values.sum()) / 3.0f;
        return sp;
    }
}

bool ReflectanceModelUtility::insertAngle0(const ReflectanceModel&  model,
                                           Brdf*                    brdf,
                                           int                      numAngles,
                                           DataType                 dataType)
{
    SampleSet* ss = brdf->getSampleSet();

    float newAngle = 0.0f;
    float maxDiff = 0.0f;
    Arrayf& angles = ss->getAngles0();

    Arrayf& angles1 = ss->getAngles1();
    Arrayf& angles2 = ss->getAngles2();
    Arrayf& angles3 = ss->getAngles3();
    for (int index = 0; index < angles.size() - 1; ++index) {
        for (int i1 = 0; i1 < angles1.size(); ++i1) {
        for (int i2 = 0; i2 < angles2.size(); ++i2) {
        for (int i3 = 0; i3 < angles3.size(); ++i3) {
            // Avoid insertion considering total internal reflection.
            auto specBrdf = dynamic_cast<SpecularCoordinatesBrdf*>(brdf);
            if (specBrdf &&
                specBrdf->getNumSpecularOffsets() != 0) {
                if (angles[index] + specBrdf->getSpecularOffset(index) > decrease(PI_2_F)) {
                    continue;
                }
            }

            int nextIndex = index + 1;

            if (hasDownwardDir(*brdf, index,     i1, i2, i3) &&
                hasDownwardDir(*brdf, nextIndex, i1, i2, i3)) break;

            float angle     = angles[index];
            float nextAngle = angles[nextIndex];
            float midAngle  = (angle + nextAngle) * 0.5f;

            float interval = nextAngle - angle;

            Spectrum sp     = getSpectrum(model, *brdf, dataType, angle,     angles1[i1], angles2[i2], angles3[i3]);
            Spectrum nextSp = getSpectrum(model, *brdf, dataType, nextAngle, angles1[i1], angles2[i2], angles3[i3]);
            Spectrum midSp  = getSpectrum(model, *brdf, dataType, midAngle,  angles1[i1], angles2[i2], angles3[i3]);

            Spectrum interpolatedSp = (sp + nextSp) * 0.5f;
            Spectrum diffSp = (interpolatedSp - midSp).abs() * interval;

            if (maxDiff < diffSp.maxCoeff()) {
                maxDiff = diffSp.maxCoeff();
                newAngle = midAngle;
            }
        }}}
    }

    if (maxDiff < EPSILON_F) {
        return true;
    }
    else {
        array_util::appendElement(&angles, newAngle);
        std::sort(angles.data(), angles.data() + angles.size());

        if (angles.size() == numAngles) {
            return true;
        }
    }

    return false;
}

bool ReflectanceModelUtility::insertAngle1(const ReflectanceModel&  model,
                                           Brdf*                    brdf,
                                           int                      numAngles,
                                           DataType                 dataType)
{
    SampleSet* ss = brdf->getSampleSet();

    float newAngle = 0.0f;
    float maxDiff = 0.0f;
    Arrayf& angles = ss->getAngles1();

    Arrayf& angles0 = ss->getAngles0();
    Arrayf& angles2 = ss->getAngles2();
    Arrayf& angles3 = ss->getAngles3();
    for (int index = 0; index < angles.size() - 1; ++index) {
        for (int i0 = 0; i0 < angles0.size(); ++i0) {
        for (int i2 = 0; i2 < angles2.size(); ++i2) {
        for (int i3 = 0; i3 < angles3.size(); ++i3) {
            // Avoid insertion considering total internal reflection.
            auto specBrdf = dynamic_cast<SpecularCoordinatesBrdf*>(brdf);
            if (specBrdf &&
                specBrdf->getNumSpecularOffsets() != 0) {
                if (angles0[i0] + specBrdf->getSpecularOffset(i0) > decrease(PI_2_F)) {
                    continue;
                }
            }

            int nextIndex = index + 1;

            if (hasDownwardDir(*brdf, i0, index,     i2, i3) &&
                hasDownwardDir(*brdf, i0, nextIndex, i2, i3)) break;

            float angle     = angles[index];
            float nextAngle = angles[nextIndex];
            float midAngle  = (angle + nextAngle) * 0.5f;

            float interval = nextAngle - angle;

            // Avoid the precision error of 32-bit float.
            if (interval < 0.0001f) continue;

            Spectrum sp     = getSpectrum(model, *brdf, dataType, angles0[i0], angle,     angles2[i2], angles3[i3]);
            Spectrum nextSp = getSpectrum(model, *brdf, dataType, angles0[i0], nextAngle, angles2[i2], angles3[i3]);
            Spectrum midSp  = getSpectrum(model, *brdf, dataType, angles0[i0], midAngle,  angles2[i2], angles3[i3]);

            Spectrum interpolatedSp = (sp + nextSp) * 0.5f;
            Spectrum diffSp = (interpolatedSp - midSp).abs() * interval;

            if (maxDiff < diffSp.maxCoeff()) {
                maxDiff = diffSp.maxCoeff();
                newAngle = midAngle;
            }
        }}}
    }

    if (maxDiff < EPSILON_F) {
        return true;
    }
    else {
        array_util::appendElement(&angles, newAngle);
        std::sort(angles.data(), angles.data() + angles.size());

        if (angles.size() == numAngles) {
            return true;
        }
    }

    return false;
}

bool ReflectanceModelUtility::insertAngle2(const ReflectanceModel&  model,
                                           Brdf*                    brdf,
                                           int                      numAngles,
                                           DataType                 dataType)
{
    SampleSet* ss = brdf->getSampleSet();

    float newAngle = 0.0f;
    float maxDiff = 0.0f;
    Arrayf& angles = ss->getAngles2();

    Arrayf& angles0 = ss->getAngles0();
    Arrayf& angles1 = ss->getAngles1();
    Arrayf& angles3 = ss->getAngles3();
    for (int index = 0; index < angles.size() - 1; ++index) {
        for (int i0 = 0; i0 < angles0.size(); ++i0) {
        for (int i1 = 0; i1 < angles1.size(); ++i1) {
        for (int i3 = 0; i3 < angles3.size(); ++i3) {
            // Avoid insertion considering total internal reflection.
            auto specBrdf = dynamic_cast<SpecularCoordinatesBrdf*>(brdf);
            if (specBrdf &&
                specBrdf->getNumSpecularOffsets() != 0) {
                if (angles0[i0] + specBrdf->getSpecularOffset(i0) > decrease(PI_2_F)) {
                    continue;
                }
            }

            int nextIndex = index + 1;

            if (hasDownwardDir(*brdf, i0, i1, index,     i3) &&
                hasDownwardDir(*brdf, i0, i1, nextIndex, i3)) break;

            float angle     = angles[index];
            float nextAngle = angles[nextIndex];
            float midAngle  = (angle + nextAngle) * 0.5f;

            float interval = nextAngle - angle;

            Spectrum sp     = getSpectrum(model, *brdf, dataType, angles0[i0], angles1[i1], angle,     angles3[i3]);
            Spectrum nextSp = getSpectrum(model, *brdf, dataType, angles0[i0], angles1[i1], nextAngle, angles3[i3]);
            Spectrum midSp  = getSpectrum(model, *brdf, dataType, angles0[i0], angles1[i1], midAngle,  angles3[i3]);

            Spectrum interpolatedSp = (sp + nextSp) * 0.5f;
            Spectrum diffSp = (interpolatedSp - midSp).abs() * interval;

            if (maxDiff < diffSp.maxCoeff()) {
                maxDiff = diffSp.maxCoeff();
                newAngle = midAngle;
            }
        }}}
    }

    if (maxDiff < EPSILON_F) {
        return true;
    }
    else {
        array_util::appendElement(&angles, newAngle);
        std::sort(angles.data(), angles.data() + angles.size());

        if (angles.size() == numAngles) {
            return true;
        }
    }

    return false;
}

bool ReflectanceModelUtility::insertAngle3(const ReflectanceModel&  model,
                                           Brdf*                    brdf,
                                           int                      numAngles,
                                           DataType                 dataType)
{
    SampleSet* ss = brdf->getSampleSet();

    float newAngle = 0.0f;
    float maxDiff = 0.0f;
    Arrayf& angles = ss->getAngles3();

    Arrayf& angles0 = ss->getAngles0();
    Arrayf& angles1 = ss->getAngles1();
    Arrayf& angles2 = ss->getAngles2();
    for (int index = 0; index < angles.size() - 1; ++index) {
        for (int i0 = 0; i0 < angles0.size(); ++i0) {
        for (int i1 = 0; i1 < angles1.size(); ++i1) {
        for (int i2 = 0; i2 < angles2.size(); ++i2) {
            // Avoid insertion considering total internal reflection.
            auto specBrdf = dynamic_cast<SpecularCoordinatesBrdf*>(brdf);
            if (specBrdf &&
                specBrdf->getNumSpecularOffsets() != 0) {
                if (angles0[i0] + specBrdf->getSpecularOffset(i0) > decrease(PI_2_F)) {
                    continue;
                }
            }

            int nextIndex = index + 1;

            if (hasDownwardDir(*brdf, i0, i1, i2, index) &&
                hasDownwardDir(*brdf, i0, i1, i2, nextIndex)) break;

            float angle     = angles[index];
            float nextAngle = angles[nextIndex];
            float midAngle  = (angle + nextAngle) * 0.5f;

            float interval = nextAngle - angle;

            // Avoid the precision error of 32-bit float.
            if (interval < 0.0001f) continue;

            Spectrum sp     = getSpectrum(model, *brdf, dataType, angles0[i0], angles1[i1], angles2[i2], angle);
            Spectrum nextSp = getSpectrum(model, *brdf, dataType, angles0[i0], angles1[i1], angles2[i2], nextAngle);
            Spectrum midSp  = getSpectrum(model, *brdf, dataType, angles0[i0], angles1[i1], angles2[i2], midAngle);

            Spectrum interpolatedSp = (sp + nextSp) * 0.5f;
            Spectrum diffSp = (interpolatedSp - midSp).abs() * interval;

            if (maxDiff < diffSp.maxCoeff()) {
                maxDiff = diffSp.maxCoeff();
                newAngle = midAngle;
            }
        }}}
    }

    if (maxDiff < EPSILON_F) {
        return true;
    }
    else {
        array_util::appendElement(&angles, newAngle);
        std::sort(angles.data(), angles.data() + angles.size());

        if (angles.size() == numAngles) {
            return true;
        }
    }

    return false;
}
