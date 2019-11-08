// =================================================================== //
// Copyright (C) 2018-2019 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#include <libbsdf/Brdf/Analyzer.h>

#include <libbsdf/Common/SolidAngle.h>
#include <libbsdf/ReflectanceModel/Fresnel.h>

using namespace lb;

// Private functions.
namespace {

// Compute the reflectance of a rectangle consisting of four outgoing directions.
Arrayd computeReflectanceOfRectangle(const SphericalCoordinatesBrdf& brdf,
                                     int inThIndex, int inPhIndex, int outThIndex, int outPhIndex)
{
    Spectrum sp = brdf.getSpectrum(inThIndex, inPhIndex, outThIndex,     outPhIndex)
                + brdf.getSpectrum(inThIndex, inPhIndex, outThIndex + 1, outPhIndex)
                + brdf.getSpectrum(inThIndex, inPhIndex, outThIndex,     outPhIndex + 1)
                + brdf.getSpectrum(inThIndex, inPhIndex, outThIndex + 1, outPhIndex + 1);
    sp /= 4.0f;

    double theta     = brdf.getOutTheta(outThIndex);
    double nextTheta = brdf.getOutTheta(outThIndex + 1);
    double phi       = brdf.getOutPhi(outPhIndex);
    double nextPhi   = brdf.getOutPhi(outPhIndex + 1);

    using std::cos;

    double midCosTheta = (cos(theta) + cos(nextTheta)) * 0.5;
    double solidAngle = SolidAngle::fromRectangle(theta, nextTheta, phi, nextPhi);
    return sp.cast<Arrayd::Scalar>() * midCosTheta * solidAngle;
}

}

Spectrum lb::computeReflectance(const SphericalCoordinatesBrdf& brdf, int inThIndex, int inPhIndex)
{
    Arrayd sumSpectrum;
    sumSpectrum.resize(brdf.getSampleSet()->getNumWavelengths());
    sumSpectrum.setZero();

    for (int thIndex = 0; thIndex < brdf.getNumOutTheta() - 1; ++thIndex) {
    for (int phIndex = 0; phIndex < brdf.getNumOutPhi()   - 1; ++phIndex) {
        sumSpectrum += computeReflectanceOfRectangle(brdf, inThIndex, inPhIndex, thIndex, phIndex);
    }}

    return sumSpectrum.cast<Spectrum::Scalar>();
}

Spectrum lb::computeReflectance(const SpecularCoordinatesBrdf& brdf, int inThIndex, int inPhIndex)
{
    Arrayd sumSpectrum;
    sumSpectrum.resize(brdf.getSampleSet()->getNumWavelengths());
    sumSpectrum.setZero();

    float inTheta = brdf.getInTheta(inThIndex);
    float inPhi   = brdf.getInPhi(inPhIndex);

    // An incoming polar angle of zero is offset to validate an incoming azimuthal angle.
    float offsetInTheta = std::max(inTheta, EPSILON_F);
    Vec3 inDir = SphericalCoordinateSystem::toXyz(offsetInTheta, inPhi);

    for (int thIndex = 0; thIndex < brdf.getNumSpecTheta() - 1; ++thIndex) {
    for (int phIndex = 0; phIndex < brdf.getNumSpecPhi()   - 1; ++phIndex) {
        Vec3 outDir0 = brdf.getOutDirection(inThIndex, inPhIndex, thIndex,     phIndex);
        Vec3 outDir1 = brdf.getOutDirection(inThIndex, inPhIndex, thIndex,     phIndex + 1);
        Vec3 outDir2 = brdf.getOutDirection(inThIndex, inPhIndex, thIndex + 1, phIndex + 1);
        Vec3 outDir3 = brdf.getOutDirection(inThIndex, inPhIndex, thIndex + 1, phIndex);

        Vec3 centroid;
        double solidAngle = SolidAngle::fromRectangleOnHemisphere(outDir0, outDir1, outDir2, outDir3, &centroid);

        if (solidAngle <= 0.0) continue;

        Spectrum sp = brdf.getSpectrum(inDir, centroid);
        sumSpectrum += sp.cast<Arrayd::Scalar>() * centroid.z() * solidAngle;
    }}

    return sumSpectrum.cast<Spectrum::Scalar>();
}

SampleSet2D* lb::computeReflectances(const SpecularCoordinatesBrdf& brdf)
{
    const SampleSet* ss = brdf.getSampleSet();

    SampleSet2D* reflectances = new SampleSet2D(brdf.getNumInTheta(),
                                                brdf.getNumInPhi(),
                                                ss->getColorModel(),
                                                ss->getNumWavelengths());
    reflectances->getThetaArray() = ss->getAngles0();
    reflectances->getPhiArray() = ss->getAngles1();
    reflectances->getWavelengths() = ss->getWavelengths();

    Spectrum sp;
    int inPhIndex;
    #pragma omp parallel for private(sp, inPhIndex) schedule(dynamic)
    for (int inThIndex = 0; inThIndex < brdf.getNumInTheta(); ++inThIndex) {
    for (    inPhIndex = 0; inPhIndex < brdf.getNumInPhi();   ++inPhIndex) {
        sp = computeReflectance(brdf, inThIndex, inPhIndex);
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

    if (!hasSameColor(*ss, *standardSs)) {
        lbError << "[lb::computeSpecularReflectances] Color models or wavelengths do not match.";
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

    if (!hasSameColor(*ss, *standardSs)) {
        lbError << "[lb::computeSpecularReflectances] Color models or wavelengths do not match.";
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

        // Find the maximum spectrum around the specular direction.
        for (int spThIndex = 0; spThIndex < brdf.getNumSpecTheta(); ++spThIndex) {
            if (brdf.getSpecTheta(spThIndex) > maxSpecularTheta) {
                continue;
            }

            for (int spPhIndex = 0; spPhIndex < brdf.getNumSpecPhi(); ++spPhIndex) {
                const Spectrum& sp = brdf.getSpectrum(thIndex, phIndex, spThIndex, spPhIndex);

                if (brdfSp.sum() < sp.sum()) {
                    brdfSp = sp;
                }
            }
        }

        Spectrum standardBrdfSp = standardBrdf.getSpectrum(inDir, specularDir);

        float standardRef;
        if (ior == 1.0f) {
            standardRef = 1.0f;
        }
        else {
            standardRef = fresnel(ss2->getTheta(thIndex), ior);
        }

        Spectrum refSp = brdfSp / standardBrdfSp.cwiseMax(EPSILON_F) * standardRef;
        ss2->setSpectrum(thIndex, phIndex, refSp);
    }}

    return ss2;
}

Spectrum lb::findDiffuseThresholds(const Brdf&      brdf,
                                   const double&    maxTheta)
{
    const SampleSet* ss = brdf.getSampleSet();

    Spectrum thresholds(ss->getNumWavelengths());
    thresholds.fill(std::numeric_limits<Spectrum::Scalar>::max());

    for (int i0 = 0; i0 < ss->getNumAngles0(); ++i0) {
    for (int i1 = 0; i1 < ss->getNumAngles1(); ++i1) {
    for (int i2 = 0; i2 < ss->getNumAngles2(); ++i2) {
    for (int i3 = 0; i3 < ss->getNumAngles3(); ++i3) {
        Vec3 inDir, outDir;
        brdf.getInOutDirection(i0, i1, i2, i3, &inDir, &outDir);

        double inTheta = std::acos(inDir[2]);
        double outTheta = std::acos(outDir[2]);

        if (inTheta <= maxTheta && outTheta <= maxTheta) {
            thresholds = thresholds.cwiseMin(ss->getSpectrum(i0, i1, i2, i3));
        }
    }}}}

    thresholds = thresholds.cwiseMax(0.0);

    return thresholds;
}

bool lb::isInDirDependentCoordinateSystem(const Brdf& brdf)
{
    if (dynamic_cast<const SphericalCoordinatesBrdf*>(&brdf) ||
        dynamic_cast<const SpecularCoordinatesBrdf*>(&brdf)) {
        return true;
    }
    else {
        return false;
    }
}
