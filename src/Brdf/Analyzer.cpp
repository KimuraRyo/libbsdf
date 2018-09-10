// =================================================================== //
// Copyright (C) 2018 Kimura Ryo                                       //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#include <libbsdf/Brdf/Analyzer.h>

#include <iostream>

#include <libbsdf/Brdf/Brdf.h>
#include <libbsdf/Brdf/Integrator.h>
#include <libbsdf/Brdf/SampleSet2D.h>
#include <libbsdf/Brdf/SpecularCoordinatesBrdf.h>
#include <libbsdf/Brdf/SphericalCoordinatesBrdf.h>

#include <libbsdf/Common/PoissonDiskDistributionOnSphere.h>

#include <libbsdf/ReflectanceModel/Fresnel.h>

using namespace lb;

SampleSet2D* lb::computeReflectances(const SpecularCoordinatesBrdf& brdf)
{
    const SampleSet* ss = brdf.getSampleSet();

    Integrator integrator(PoissonDiskDistributionOnSphere::NUM_SAMPLES_ON_HEMISPHERE, true);

    SampleSet2D* reflectances = new SampleSet2D(brdf.getNumInTheta(),
                                                brdf.getNumInPhi(),
                                                ss->getColorModel(),
                                                ss->getNumWavelengths());
    reflectances->getThetaArray() = ss->getAngles0();
    reflectances->getPhiArray() = ss->getAngles1();
    reflectances->getWavelengths() = ss->getWavelengths();

    for (int inThIndex = 0; inThIndex < brdf.getNumInTheta(); ++inThIndex) {
    for (int inPhIndex = 0; inPhIndex < brdf.getNumInPhi();   ++inPhIndex) {
        Vec3 inDir = SphericalCoordinateSystem::toXyz(brdf.getInTheta(inThIndex),
                                                      brdf.getInPhi(inPhIndex));
        Spectrum sp = integrator.computeReflectance(brdf, inDir);
        reflectances->setSpectrum(inThIndex, inPhIndex, sp);
    }}

    return reflectances;
}

SampleSet2D* lb::computeSpecularReflectances(const Brdf&    brdf,
                                             const Brdf&    standardBrdf,
                                             float          ior)
{
    const SampleSet* ss = brdf.getSampleSet();
    const SampleSet* standardSs = standardBrdf.getSampleSet();

    if (ss->getNumWavelengths() != standardSs->getNumWavelengths() ||
        !ss->getWavelengths().isApprox(standardSs->getWavelengths())) {
        std::cerr
            << "[lb::computeSpecularReflectances] Wavelengths do not match."
            << std::endl;
        return 0;
    }

    SampleSet2D* ss2 = new SampleSet2D(ss->getNumAngles0(),
                                       ss->getNumAngles1(),
                                       ss->getColorModel(),
                                       ss->getNumWavelengths());
    ss2->getThetaArray() = ss->getAngles0();
    ss2->getPhiArray() = ss->getAngles1();
    ss2->getWavelengths() = ss->getWavelengths();

    for (int thIndex = 0; thIndex < ss2->getNumTheta(); ++thIndex) {
    for (int phIndex = 0; phIndex < ss2->getNumPhi();   ++phIndex) {
        Vec3 inDir = ss2->getDirection(thIndex, phIndex);
        Vec3 specularDir = reflect(inDir, Vec3(0.0, 0.0, 1.0));

        Spectrum brdfSp = brdf.getSpectrum(inDir, specularDir);
        Spectrum standardBrdfSp = standardBrdf.getSpectrum(inDir, specularDir);

        float standardRef;
        if (ior == 1.0f) {
            standardRef = 1.0f;
        }
        else {
            standardRef = fresnel(ss2->getTheta(thIndex), ior);
        }

        Spectrum refSp = brdfSp / standardBrdfSp * standardRef;
        ss2->setSpectrum(thIndex, phIndex, refSp);
    }}

    return ss2;
}

SampleSet2D* lb::computeSpecularReflectances(const SpecularCoordinatesBrdf& brdf,
                                             const Brdf&                    standardBrdf,
                                             float                          ior,
                                             float                          maxSpecularTheta)
{
    const SampleSet* ss = brdf.getSampleSet();
    const SampleSet* standardSs = standardBrdf.getSampleSet();

    if (ss->getNumWavelengths() != standardSs->getNumWavelengths() ||
        !ss->getWavelengths().isApprox(standardSs->getWavelengths())) {
        std::cerr
            << "[lb::computeSpecularReflectances] Wavelengths do not match."
            << std::endl;
        return 0;
    }

    SampleSet2D* ss2 = new SampleSet2D(ss->getNumAngles0(),
                                       ss->getNumAngles1(),
                                       ss->getColorModel(),
                                       ss->getNumWavelengths());
    ss2->getThetaArray()    = ss->getAngles0();
    ss2->getPhiArray()      = ss->getAngles1();
    ss2->getWavelengths()   = ss->getWavelengths();

    for (int thIndex = 0; thIndex < ss2->getNumTheta(); ++thIndex) {
    for (int phIndex = 0; phIndex < ss2->getNumPhi();   ++phIndex) {
        Vec3 inDir = ss2->getDirection(thIndex, phIndex);
        Vec3 specularDir = reflect(inDir, Vec3(0.0, 0.0, 1.0));

        Spectrum brdfSp = brdf.getSpectrum(inDir, specularDir);

        for (int spThIndex = 0; spThIndex < brdf.getNumSpecTheta(); ++spThIndex) {
        for (int spPhIndex = 0; spPhIndex < brdf.getNumSpecPhi();   ++spPhIndex) {
            const Spectrum& sp = brdf.getSpectrum(thIndex, phIndex, spThIndex, spPhIndex);

            if (brdfSp.sum() < sp.sum()) {
                brdfSp = sp;
            }
        }}

        Spectrum standardBrdfSp = standardBrdf.getSpectrum(inDir, specularDir);

        float standardRef;
        if (ior == 1.0f) {
            standardRef = 1.0f;
        }
        else {
            standardRef = fresnel(ss2->getTheta(thIndex), ior);
        }

        Spectrum refSp = brdfSp / standardBrdfSp * standardRef;
        ss2->setSpectrum(thIndex, phIndex, refSp);
    }}

    return ss2;
}


Spectrum lb::findDiffuseThresholds(const lb::Brdf&  brdf,
                                   float            maxTheta)
{
    const lb::SampleSet* ss = brdf.getSampleSet();

    lb::Spectrum thresholds(ss->getNumWavelengths());
    thresholds.fill(std::numeric_limits<lb::Spectrum::Scalar>::max());

    for (int i0 = 0; i0 < ss->getNumAngles0(); ++i0) {
    for (int i1 = 0; i1 < ss->getNumAngles1(); ++i1) {
    for (int i2 = 0; i2 < ss->getNumAngles2(); ++i2) {
    for (int i3 = 0; i3 < ss->getNumAngles3(); ++i3) {
        lb::Vec3 inDir, outDir;
        brdf.getInOutDirection(i0, i1, i2, i3, &inDir, &outDir);

        float inTheta = std::acos(inDir[2]);
        float outTheta = std::acos(outDir[2]);

        if (inTheta <= maxTheta && outTheta <= maxTheta) {
            thresholds = thresholds.cwiseMin(ss->getSpectrum(i0, i1, i2, i3));
        }
    }}}}

    thresholds = thresholds.cwiseMax(0.0);

    return thresholds;
}
