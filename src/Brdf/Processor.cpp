// =================================================================== //
// Copyright (C) 2014-2015 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#include <libbsdf/Brdf/Processor.h>

#include <iostream>

#include <libbsdf/Brdf/Integrator.h>
#include <libbsdf/Brdf/RandomSampleSet.h>
#include <libbsdf/Brdf/SampleSet2D.h>

#include <libbsdf/Common/PoissonDiskDistributionOnSphere.h>
#include <libbsdf/Common/SpectrumUtility.h>
#include <libbsdf/Common/SphericalCoordinateSystem.h>

using namespace lb;

bool lb::initizlizeSpectra(const Brdf& baseBrdf, Brdf* brdf)
{
    const SampleSet* baseSs = baseBrdf.getSampleSet();
    SampleSet* ss = brdf->getSampleSet();

    bool same = true;

    if (baseSs->getColorModel() != ss->getColorModel()) {
        same = false;
        std::cerr
            << "[initizlizeSpectra] Color models don't match: "
            << baseSs->getColorModel() << ", " << ss->getColorModel()
            << std::endl;
    }

    if (!baseSs->getWavelengths().isApprox(ss->getWavelengths())) {
        same = false;
        std::cerr
            << "[initizlizeSpectra] Wavelengths don't match: "
            << baseSs->getWavelengths() << ", " << ss->getWavelengths()
            << std::endl;
    }

    if (!same) return false;

    for (int i0 = 0; i0 < ss->getNumAngles0(); ++i0) {
    for (int i1 = 0; i1 < ss->getNumAngles1(); ++i1) {
    for (int i2 = 0; i2 < ss->getNumAngles2(); ++i2) {
    for (int i3 = 0; i3 < ss->getNumAngles3(); ++i3) {
        Vec3 inDir, outDir;
        brdf->getInOutDirection(i0, i1, i2, i3, &inDir, &outDir);

        ss->setSpectrum(i0, i1, i2, i3, baseBrdf.getSpectrum(inDir, outDir));
    }}}}

    return true;
}

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
        if (cosOutTheta <= 0.0f && i2 > 0) {
            // Assume i2 is the index of the polar angle related to outgoing directions.
            brdf->getInOutDirection(i0, i1, i2 - 1, i3, &inDir, &outDir);
            sp = ss->getSpectrum(i0, i1, i2 - 1, i3);
            cosOutTheta = outDir.dot(Vec3(0.0, 0.0, 1.0));
        }

        sp /= cosOutTheta;
    }}}}
}

SphericalCoordinatesBrdf* lb::fillSymmetricBrdf(SphericalCoordinatesBrdf* brdf)
{
    RandomSampleSet::AngleList filledAngles;

    for (int i = 0; i < brdf->getNumOutPhi(); ++i) {
        float outPhi = brdf->getOutPhi(i);
        bool angleOmitted = (outPhi != 0.0f &&
                             !isEqual(outPhi, PI_F) &&
                             !isEqual(outPhi, 2.0f * PI_F));
        if (angleOmitted) {
            filledAngles.push_back(SphericalCoordinateSystem::MAX_ANGLE3 - outPhi);
        }
    }

    const SampleSet* ss = brdf->getSampleSet();

    SphericalCoordinatesBrdf* filledBrdf = new SphericalCoordinatesBrdf(brdf->getNumInTheta(),
                                                                        brdf->getNumInPhi(),
                                                                        brdf->getNumOutTheta(),
                                                                        brdf->getNumOutPhi() + filledAngles.size(),
                                                                        ss->getColorModel(),
                                                                        ss->getNumWavelengths());
    SampleSet* filledSs = filledBrdf->getSampleSet();

    // Set angles.
    filledSs->getAngles0() = ss->getAngles0();
    filledSs->getAngles1() = ss->getAngles1();
    filledSs->getAngles2() = ss->getAngles2();
    for (int i = 0; i < filledBrdf->getNumOutPhi(); ++i) {
        if (i < brdf->getNumOutPhi()) {
            filledBrdf->setOutPhi(i, brdf->getOutPhi(i));
        }
        else {
            filledBrdf->setOutPhi(i, filledAngles.at(i - brdf->getNumOutPhi()));
        }
    }
    Arrayf& outPhiAngles = filledSs->getAngles3();
    std::sort(outPhiAngles.data(), outPhiAngles.data() + outPhiAngles.size());

    // Set wavelengths.
    for (int i = 0; i < filledSs->getNumWavelengths(); ++i) {
        float wl = ss->getWavelength(i);
        filledSs->setWavelength(i, wl);
    }

    for (int inThIndex  = 0; inThIndex  < filledBrdf->getNumInTheta();  ++inThIndex)  {
    for (int inPhIndex  = 0; inPhIndex  < filledBrdf->getNumInPhi();    ++inPhIndex)  {
    for (int outThIndex = 0; outThIndex < filledBrdf->getNumOutTheta(); ++outThIndex) {
    for (int outPhIndex = 0; outPhIndex < filledBrdf->getNumOutPhi();   ++outPhIndex) {
        float outPhi = filledBrdf->getOutPhi(outPhIndex);

        // Find the corresponding index.
        int origIndex;
        for (origIndex = 0; origIndex < brdf->getNumOutPhi(); ++origIndex) {
            float origOutPhi = brdf->getOutPhi(origIndex);
            bool outPhiEqual = (origOutPhi == outPhi ||
                                isEqual(origOutPhi, SphericalCoordinateSystem::MAX_ANGLE3 - outPhi));
            if (outPhiEqual) break;
        }

        Spectrum& sp = brdf->getSpectrum(inThIndex, inPhIndex, outThIndex, origIndex);
        filledBrdf->setSpectrum(inThIndex, inPhIndex, outThIndex, outPhIndex, sp);
    }}}}

    return filledBrdf;
}

SphericalCoordinatesBrdf* lb::rotateOutPhi(const SphericalCoordinatesBrdf&  brdf,
                                           float                            rotationAngle)
{
    assert(rotationAngle > -2.0f * PI_F && rotationAngle < 2.0f * PI_F);

    if (rotationAngle < 0.0f) {
        rotationAngle += 2.0f * PI_F;
    }

    SphericalCoordinatesBrdf* rotatedBrdf = new SphericalCoordinatesBrdf(brdf);
    SampleSet* ss = rotatedBrdf->getSampleSet();

    ss->updateAngleAttributes();
    if (!ss->isEqualIntervalAngles3()) {
        for (int i = 0; i < rotatedBrdf->getNumOutPhi(); ++i) {
            float outPhi = rotatedBrdf->getOutPhi(i) + rotationAngle;
            if (outPhi > 2.0f * PI_F) {
                outPhi -= 2.0f * PI_F;
            }

            rotatedBrdf->setOutPhi(i, outPhi);
        }

        Arrayf& outPhiAngles = ss->getAngles3();
        std::sort(outPhiAngles.data(), outPhiAngles.data() + outPhiAngles.size());
    }

    for (int inThIndex  = 0; inThIndex  < rotatedBrdf->getNumInTheta();  ++inThIndex)  {
    for (int inPhIndex  = 0; inPhIndex  < rotatedBrdf->getNumInPhi();    ++inPhIndex)  {
    for (int outThIndex = 0; outThIndex < rotatedBrdf->getNumOutTheta(); ++outThIndex) {
    for (int outPhIndex = 0; outPhIndex < rotatedBrdf->getNumOutPhi();   ++outPhIndex) {
        float inTheta  = rotatedBrdf->getInTheta(inThIndex);
        float inPhi    = rotatedBrdf->getInPhi(inPhIndex);
        float outTheta = rotatedBrdf->getOutTheta(outThIndex);
        float outPhi   = rotatedBrdf->getOutPhi(outPhIndex) - rotationAngle;

        if (outPhi < 0.0f) {
            outPhi += 2.0f * PI_F;
        }

        Spectrum sp = brdf.getSpectrum(inTheta, inPhi, outTheta, outPhi);
        rotatedBrdf->setSpectrum(inThIndex, inPhIndex, outThIndex, outPhIndex, sp);
    }}}}

    return rotatedBrdf;
}

void lb::fixEnergyConservation(SpecularCoordinatesBrdf* brdf)
{
    SampleSet* ss = brdf->getSampleSet();

    int numInTheta = brdf->getNumInTheta();
    int numInPhi = brdf->getNumInPhi();

    SampleSet2D reflectances(numInTheta, numInPhi,
                             ss->getColorModel(), ss->getNumWavelengths());
    reflectances.getThetaArray()  = ss->getAngles0();
    reflectances.getPhiArray()    = ss->getAngles1();
    reflectances.getWavelengths() = ss->getWavelengths();

    Integrator integrator(PoissonDiskDistributionOnSphere::NUM_SAMPLES_ON_HEMISPHERE, true);

    for (int inThIndex = 0; inThIndex < numInTheta; ++inThIndex) {
    for (int inPhIndex = 0; inPhIndex < numInPhi;   ++inPhIndex) {
        Vec3 inDir = SphericalCoordinateSystem::toXyz(brdf->getInTheta(inThIndex),
                                                      brdf->getInPhi(inPhIndex));
        Spectrum sp = integrator.computeReflectance(*brdf, inDir);
        reflectances.setSpectrum(inThIndex, inPhIndex, sp);

        // Fix samples to conserve energy.
        float maxReflectance = sp.maxCoeff();
        if (maxReflectance > 1.0f) {
            for (int i2 = 0; i2 < ss->getNumAngles2(); ++i2) {
            for (int i3 = 0; i3 < ss->getNumAngles3(); ++i3) {
                Spectrum& fixedSp = ss->getSpectrum(inThIndex, inPhIndex, i2, i3);
                const float coeff = 0.999546f; // Reflectance of Lambertian using lb::Integrator.
                fixedSp /= maxReflectance / coeff;
            }}
        }
    }}
}

void lb::copySpectraFromPhiOfZeroTo2PI(Brdf* brdf)
{
    SampleSet* ss = brdf->getSampleSet();

    if (ss->getNumAngles1() >= 2 &&
        ss->getAngle1(0) == 0.0f &&
        ss->getAngle1(ss->getNumAngles1() - 1) >= SphericalCoordinateSystem::MAX_ANGLE1) {
        for (int i0 = 0; i0 < ss->getNumAngles0(); ++i0) {
        for (int i2 = 0; i2 < ss->getNumAngles2(); ++i2) {
        for (int i3 = 0; i3 < ss->getNumAngles3(); ++i3) {
            const Spectrum& sp = ss->getSpectrum(i0, 0, i2, i3);
            ss->setSpectrum(i0, ss->getNumAngles1() - 1, i2, i3, sp);
        }}}
    }

    if (ss->getNumAngles3() >= 2 &&
        ss->getAngle3(0) == 0.0f &&
        ss->getAngle3(ss->getNumAngles3() - 1) >= SphericalCoordinateSystem::MAX_ANGLE3) {
        for (int i0 = 0; i0 < ss->getNumAngles0(); ++i0) {
        for (int i1 = 0; i1 < ss->getNumAngles1(); ++i1) {
        for (int i2 = 0; i2 < ss->getNumAngles2(); ++i2) {
            const Spectrum& sp = ss->getSpectrum(i0, i1, i2, 0);
            ss->setSpectrum(i0, i1, i2, ss->getNumAngles3() - 1, sp);
        }}}
    }
}

void lb::convertFromXyzToSrgb(SampleSet* samples)
{
    ColorModel cm = samples->getColorModel();
    if (cm != XYZ_MODEL) {
        std::cerr << "[convertFromXyzToSrgb] Not CIE-XYZ model: " << cm << std::endl;
        return;
    }

    for (int i0 = 0; i0 < samples->getNumAngles0(); ++i0) {
    for (int i1 = 0; i1 < samples->getNumAngles1(); ++i1) {
    for (int i2 = 0; i2 < samples->getNumAngles2(); ++i2) {
    for (int i3 = 0; i3 < samples->getNumAngles3(); ++i3) {
        const Spectrum& xyz = samples->getSpectrum(i0, i1, i2, i3);
        Spectrum rgb = SpectrumUtility::xyzToSrgb(xyz);
        samples->setSpectrum(i0, i1, i2, i3, rgb);
    }}}}

    samples->setColorModel(RGB_MODEL);
}

void lb::fillSpectra(SampleSet* samples, float value)
{
    lb::SpectrumList& spectra = samples->getSpectra();
    for (auto it = spectra.begin(); it != spectra.end(); ++it) {
        it->fill(value);
    }
}

void lb::multiplySpectra(SampleSet* samples, float value)
{
    lb::SpectrumList& spectra = samples->getSpectra();
    for (auto it = spectra.begin(); it != spectra.end(); ++it) {
        *it *= value;
    }
}

void lb::clampNegativeSpectra(SampleSet* samples)
{
    lb::SpectrumList& spectra = samples->getSpectra();
    for (auto it = spectra.begin(); it != spectra.end(); ++it) {
        *it = it->cwiseMax(0.0f);
    }
}
