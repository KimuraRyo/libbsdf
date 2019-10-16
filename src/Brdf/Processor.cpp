// =================================================================== //
// Copyright (C) 2014-2019 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#include <libbsdf/Brdf/Processor.h>

#include <libbsdf/Brdf/Analyzer.h>
#include <libbsdf/Brdf/RandomSampleSet.h>

using namespace lb;

void lb::editComponents(const Brdf&         origBrdf,
                        Brdf*               brdf,
                        const Spectrum&     diffuseThresholds,
                        Spectrum::Scalar    glossyIntensity,
                        Spectrum::Scalar    glossyShininess,
                        Spectrum::Scalar    diffuseIntensity)
{
    const SampleSet* ss = brdf->getSampleSet();

    for (int i0 = 0; i0 < ss->getNumAngles0(); ++i0) {
    for (int i1 = 0; i1 < ss->getNumAngles1(); ++i1) {
    for (int i2 = 0; i2 < ss->getNumAngles2(); ++i2) {
    #pragma omp parallel for
    for (int i3 = 0; i3 < ss->getNumAngles3(); ++i3) {
        editComponents(i0, i1, i2, i3,
                       origBrdf, brdf,
                       diffuseThresholds,
                       glossyIntensity, glossyShininess, diffuseIntensity);
    }}}}

    brdf->setSourceType(EDITED_SOURCE);
}

void lb::editComponents(int                 i0,
                        int                 i1,
                        int                 i2,
                        int                 i3,
                        const Brdf&         origBrdf,
                        Brdf*               brdf,
                        const Spectrum&     diffuseThresholds,
                        Spectrum::Scalar    glossyIntensity,
                        Spectrum::Scalar    glossyShininess,
                        Spectrum::Scalar    diffuseIntensity)
{
    using std::pow;
    using std::max;

    SampleSet* ss = brdf->getSampleSet();

    Vec3 inDir, outDir;
    brdf->getInOutDirection(i0, i1, i2, i3, &inDir, &outDir);

    // Offset outgoing directions to edit a glossy component with shininess.
    if (glossyShininess != 1.0) {
        float inTh, inPh, specTh, specPh;
        SpecularCoordinateSystem::fromXyz(inDir, outDir, &inTh, &inPh, &specTh, &specPh);
        float specThWeight = specTh / SpecularCoordinateSystem::MAX_ANGLE2;
        specThWeight = pow(specThWeight, 1.0f / max(glossyShininess, EPSILON_F));
        float newSpecTh = specThWeight * SpecularCoordinateSystem::MAX_ANGLE2;

        inTh        = clamp(inTh,       SpecularCoordinateSystem::MIN_ANGLE0, SpecularCoordinateSystem::MAX_ANGLE0);
        inPh        = clamp(inPh,       SpecularCoordinateSystem::MIN_ANGLE1, SpecularCoordinateSystem::MAX_ANGLE1);
        newSpecTh   = clamp(newSpecTh,  SpecularCoordinateSystem::MIN_ANGLE2, SpecularCoordinateSystem::MAX_ANGLE2);
        specPh      = clamp(specPh,     SpecularCoordinateSystem::MIN_ANGLE3, SpecularCoordinateSystem::MAX_ANGLE3);

        SpecularCoordinateSystem::toXyz(inTh, inPh, newSpecTh, specPh, &inDir, &outDir);
    }

    Spectrum origSp = origBrdf.getSpectrum(inDir, outDir);
    Spectrum sp(origSp.size());

    // Edit a BRDF with glossy and diffuse intensity.
    for (int i = 0; i < sp.size(); ++i) {
        Spectrum::Scalar origVal = origSp[i];
        Spectrum::Scalar threshold = diffuseThresholds[i];
        if (origVal <= threshold) {
            sp[i] = origVal * diffuseIntensity;
        }
        else {
            Spectrum::Scalar glossy = origVal - threshold;
            sp[i] = glossy * glossyIntensity + threshold * diffuseIntensity;
        }
    }

    ss->setSpectrum(i0, i1, i2, i3, sp);
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
        Vec3::Scalar cosOutTheta = outDir.dot(Vec3(0.0, 0.0, 1.0));

        Spectrum& sp = ss->getSpectrum(i0, i1, i2, i3);

        // Copy the spectrum if the Z-component of the outgoing direction is zero or negative.
        if (cosOutTheta <= 0.0 && i2 > 0) {
            // Assume i2 is the index of the polar angle related to outgoing directions.
            brdf->getInOutDirection(i0, i1, i2 - 1, i3, &inDir, &outDir);
            sp = ss->getSpectrum(i0, i1, i2 - 1, i3);
            cosOutTheta = outDir.dot(Vec3(0.0, 0.0, 1.0));
        }

        sp /= static_cast<Spectrum::Scalar>(cosOutTheta);
    }}}}
}

SphericalCoordinatesBrdf* lb::fillSymmetricBrdf(SphericalCoordinatesBrdf* brdf)
{
    RandomSampleSet<SphericalCoordinateSystem>::AngleList filledAngles;

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

    int numOutPhi = brdf->getNumOutPhi() + static_cast<int>(filledAngles.size());
    SphericalCoordinatesBrdf* filledBrdf = new SphericalCoordinatesBrdf(brdf->getNumInTheta(),
                                                                        brdf->getNumInPhi(),
                                                                        brdf->getNumOutTheta(),
                                                                        numOutPhi,
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

void lb::fillSpectraAtInThetaOf0(Brdf* brdf)
{
    SampleSet* ss = brdf->getSampleSet();

    bool acceptable = (dynamic_cast<SphericalCoordinatesBrdf*>(brdf) ||
                       dynamic_cast<SpecularCoordinatesBrdf*>(brdf));

    if (!acceptable ||
        !ss->isIsotropic() ||
        ss->getAngle0(0) != 0.0f) {
        return;
    }

    for (int i2 = 0; i2 < ss->getNumAngles2(); ++i2) {
        int numSamples = 0;
        Spectrum sumSp = Spectrum::Zero(ss->getNumWavelengths());

        for (int i3 = 0; i3 < ss->getNumAngles3(); ++i3) {
            if (i3 == ss->getNumAngles3() - 1 &&
                isEqual(ss->getAngle3(0), ss->getAngle3(ss->getNumAngles3() - 1)) &&
                ss->getNumAngles3() > 1) {
                break;
            }

            sumSp += ss->getSpectrum(0, i2, i3);
            ++numSamples;
        }

        Spectrum avgSp = sumSp / static_cast<Spectrum::Scalar>(numSamples);

        for (int i3 = 0; i3 < ss->getNumAngles3(); ++i3) {
            ss->setSpectrum(0, i2, i3, avgSp);
        }
    }
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

    fixNegativeSpectra(ss);

    for (int inThIndex = 0; inThIndex < brdf->getNumInTheta(); ++inThIndex) {
    for (int inPhIndex = 0; inPhIndex < brdf->getNumInPhi();   ++inPhIndex) {
        Spectrum sp = computeReflectance(*brdf, inThIndex, inPhIndex);

        // Fix samples to conserve energy.
        float maxReflectance = sp.maxCoeff();
        if (maxReflectance > 1.0f) {
            for (int i2 = 0; i2 < ss->getNumAngles2(); ++i2) {
            for (int i3 = 0; i3 < ss->getNumAngles3(); ++i3) {
                Spectrum& fixedSp = ss->getSpectrum(inThIndex, inPhIndex, i2, i3);
                fixedSp /= maxReflectance;
            }}
        }
    }}
}

void lb::fixEnergyConservation(SpecularCoordinatesBrdf* brdf,
                               const SampleSet2D&       specularReflectances)
{
    SampleSet* ss = brdf->getSampleSet();

    for (int inThIndex = 0; inThIndex < brdf->getNumInTheta(); ++inThIndex) {
    for (int inPhIndex = 0; inPhIndex < brdf->getNumInPhi();   ++inPhIndex) {
        Spectrum sp = computeReflectance(*brdf, inThIndex, inPhIndex);

        // Fix samples to conserve energy.
        float maxReflectance = sp.maxCoeff();
        if (maxReflectance > 1.0f) {
            for (int i2 = 0; i2 < ss->getNumAngles2(); ++i2) {
            for (int i3 = 0; i3 < ss->getNumAngles3(); ++i3) {
                Spectrum& fixedSp = ss->getSpectrum(inThIndex, inPhIndex, i2, i3);
                fixedSp /= maxReflectance;
            }}
        }

        Spectrum specRefSp = specularReflectances.getSpectrum(brdf->getInTheta(inThIndex),
                                                              brdf->getInPhi(inPhIndex));
        sp += specRefSp;

        maxReflectance = 0;
        int maxIndex = 0;
        for (int i = 0; i < sp.size(); ++i) {
            if (sp[i] > maxReflectance) {
                maxReflectance = sp[i];
                maxIndex = i;
            }
        }

        // Fix samples to conserve energy with specular reflectances.
        if (maxReflectance > 1.0f) {
            for (int i2 = 0; i2 < ss->getNumAngles2(); ++i2) {
            for (int i3 = 0; i3 < ss->getNumAngles3(); ++i3) {
                Spectrum& fixedSp = ss->getSpectrum(inThIndex, inPhIndex, i2, i3);
                fixedSp *= 1.0f - specRefSp[maxIndex];
            }}
        }
    }}
}

void lb::fixEnergyConservation(SpecularCoordinatesBrdf* brdf,
                               SpecularCoordinatesBrdf* btdf)
{
    fixNegativeSpectra(brdf->getSampleSet());
    fixNegativeSpectra(btdf->getSampleSet());

    SampleSet2D* reflectances = computeReflectances(*brdf);
    SampleSet2D* transmittances = computeReflectances(*btdf);

    // Process BRDF.
    for (int inThIndex = 0; inThIndex < brdf->getNumInTheta(); ++inThIndex) {
    for (int inPhIndex = 0; inPhIndex < brdf->getNumInPhi();   ++inPhIndex) {
        Spectrum sp = reflectances->getSpectrum(inThIndex, inPhIndex)
                    + transmittances->getSpectrum(btdf->getInTheta(inThIndex), btdf->getInPhi(inPhIndex));

        // Fix samples to conserve energy.
        float maxReflectance = sp.maxCoeff();
        if (maxReflectance > 1.0f) {
            for (int spThIndex = 0; spThIndex < brdf->getNumSpecTheta(); ++spThIndex) {
            for (int spPhIndex = 0; spPhIndex < brdf->getNumSpecPhi();   ++spPhIndex) {
                Spectrum& fixedSp = brdf->getSpectrum(inThIndex, inPhIndex, spThIndex, spPhIndex);
                fixedSp /= maxReflectance;
            }}
        }
    }}

    // Process BTDF.
    for (int inThIndex = 0; inThIndex < btdf->getNumInTheta(); ++inThIndex) {
    for (int inPhIndex = 0; inPhIndex < btdf->getNumInPhi();   ++inPhIndex) {
        Spectrum sp = reflectances->getSpectrum(btdf->getInTheta(inThIndex), btdf->getInPhi(inPhIndex))
                    + transmittances->getSpectrum(inThIndex, inPhIndex);

        // Fix samples to conserve energy.
        float maxReflectance = sp.maxCoeff();
        if (maxReflectance > 1.0f) {
            for (int spThIndex = 0; spThIndex < btdf->getNumSpecTheta(); ++spThIndex) {
            for (int spPhIndex = 0; spPhIndex < btdf->getNumSpecPhi();   ++spPhIndex) {
                Spectrum& fixedSp = btdf->getSpectrum(inThIndex, inPhIndex, spThIndex, spPhIndex);
                fixedSp /= maxReflectance;
            }}
        }
    }}

    delete reflectances;
    delete transmittances;
}

void lb::fillBackSide(SpecularCoordinatesBrdf* brdf)
{
    for (int inThIndex = 0; inThIndex < brdf->getNumInTheta();   ++inThIndex) {
    for (int inPhIndex = 0; inPhIndex < brdf->getNumInPhi();     ++inPhIndex) {
    for (int spThIndex = 0; spThIndex < brdf->getNumSpecTheta(); ++spThIndex) {
        bool spPhBoundaryFound = false;
        int boundary0 = 0;

        if (brdf->getInTheta(inThIndex) > 0.0f) {
            // Search the boundary of specular azimuthal angles.
            bool upwardDirFound = false;
            for (int spPhIndex = 0; spPhIndex < brdf->getNumSpecPhi(); ++spPhIndex) {
                Vec3 inDir, outDir;
                brdf->toXyz(brdf->getInTheta(inThIndex),
                            brdf->getInPhi(inPhIndex),
                            brdf->getSpecTheta(spThIndex),
                            brdf->getSpecPhi(spPhIndex),
                            &inDir, &outDir);

                if (!isDownwardDir(outDir)) {
                    upwardDirFound = true;
                }
                else if (upwardDirFound) {
                    spPhBoundaryFound = true;
                    break;
                }

                boundary0 = spPhIndex;
            }
        }

        if (spPhBoundaryFound) {
            int boundary1 = 0;

            // Search another boundary from the opposite direction.
            bool upwardDirFound = false;
            for (int spPhIndex = brdf->getNumSpecPhi() - 1; spPhIndex >= 0; --spPhIndex) {
                Vec3 inDir, outDir;
                brdf->toXyz(brdf->getInTheta(inThIndex),
                            brdf->getInPhi(inPhIndex),
                            brdf->getSpecTheta(spThIndex),
                            brdf->getSpecPhi(spPhIndex),
                            &inDir, &outDir);

                if (!isDownwardDir(outDir)) {
                    upwardDirFound = true;
                }
                else if (upwardDirFound) {
                    break;
                }

                boundary1 = spPhIndex;
            }

            // Fill values.
            for (int spPhIndex = 0; spPhIndex < brdf->getNumSpecPhi(); ++spPhIndex) {
                Vec3 inDir, outDir;
                brdf->toXyz(brdf->getInTheta(inThIndex),
                            brdf->getInPhi(inPhIndex),
                            brdf->getSpecTheta(spThIndex),
                            brdf->getSpecPhi(spPhIndex),
                            &inDir, &outDir);

                if (!isDownwardDir(outDir)) {
                    continue;
                }

                float spPh = brdf->getSpecPhi(spPhIndex);
                float spPh0 = brdf->getSpecPhi(boundary0);
                float spPh1 = brdf->getSpecPhi(boundary1);

                int boundary;
                if ((spPh - spPh0) <= (spPh1 - spPh)) {
                    boundary = boundary0;
                }
                else {
                    boundary = boundary1;
                }

                Spectrum sp = brdf->getSpectrum(inThIndex, inPhIndex, spThIndex, boundary);
                brdf->setSpectrum(inThIndex, inPhIndex, spThIndex, spPhIndex, sp);
            }
        }
        else if (spThIndex >= 1) {
            // Copy the sample at the previous specular polar angle.
            for (int spPhIndex = 0; spPhIndex < brdf->getNumSpecPhi(); ++spPhIndex) {
                Vec3 inDir, outDir;
                brdf->toXyz(brdf->getInTheta(inThIndex),
                            brdf->getInPhi(inPhIndex),
                            brdf->getSpecTheta(spThIndex),
                            brdf->getSpecPhi(spPhIndex),
                            &inDir, &outDir);

                if (isDownwardDir(outDir)) {
                    Spectrum sp = brdf->getSpectrum(inThIndex, inPhIndex, spThIndex - 1, spPhIndex);
                    brdf->setSpectrum(inThIndex, inPhIndex, spThIndex, spPhIndex, sp);
                }
            }
        }
    }}}
}

void lb::equalizeOverlappingSamples(SpecularCoordinatesBrdf* brdf)
{
    const SampleSet* ss = brdf->getSampleSet();

    if (brdf->getNumInPhi() >= 2) {
        if (isEqual(brdf->getInTheta(0), 0.0f)) {
            const SpecularCoordinatesBrdf origBrdf = *brdf;

            // Equalize samples for incoming azimuthal angles if an incoming polar angle is 0.
            for (int inPhIndex = 0; inPhIndex < brdf->getNumInPhi();     ++inPhIndex) {
            for (int spThIndex = 0; spThIndex < brdf->getNumSpecTheta(); ++spThIndex) {
            for (int spPhIndex = 0; spPhIndex < brdf->getNumSpecPhi();   ++spPhIndex) {
                Arrayd sumSp = Arrayd::Zero(ss->getNumWavelengths());

                Vec3 outDir = brdf->getOutDirection(0, inPhIndex, spThIndex, spPhIndex);
                for (int sampledInPhIndex = 0;
                     sampledInPhIndex < brdf->getNumInPhi();
                     ++sampledInPhIndex) {

                    // An incoming polar angle of zero is offset to validate an incoming azimuthal angle.
                    const float inTheta = EPSILON_F;
                    float inPhi = brdf->getInPhi(sampledInPhIndex);

                    // Ignore an overlapping angle.
                    if (isEqual(inPhi, 2.0f * PI_F) &&
                        isEqual(brdf->getInPhi(0), 0.0f)) {
                        continue;
                    }

                    Vec3 inDir = SphericalCoordinateSystem::toXyz(inTheta, inPhi);
                    sumSp += origBrdf.getSpectrum(inDir, outDir).cast<Arrayd::Scalar>();
                }

                int numInPhi = brdf->getNumInPhi();

                // Ignore an overlapping angle.
                if (isEqual(brdf->getInPhi(0), 0.0f) &&
                    isEqual(brdf->getInPhi(brdf->getNumInPhi() - 1), 2.0f * PI_F)) {
                    --numInPhi;
                }

                Spectrum sp = sumSp.cast<Spectrum::Scalar>() / numInPhi;
                brdf->setSpectrum(0, inPhIndex, spThIndex, spPhIndex, sp);
            }}}
        }

        int minInPhiIndex = 0;
        int maxInPhiIndex = brdf->getNumInPhi() - 1;
        if (isEqual(brdf->getInPhi(minInPhiIndex), SpecularCoordinateSystem::MIN_ANGLE1) &&
            isEqual(brdf->getInPhi(maxInPhiIndex), SpecularCoordinateSystem::MAX_ANGLE1)) {
            // Equalize samples if an incoming azimuthal angle is 0 or 2PI radian.
            for (int inThIndex = 0; inThIndex < brdf->getNumInTheta();   ++inThIndex) {
            for (int spThIndex = 0; spThIndex < brdf->getNumSpecTheta(); ++spThIndex) {
            for (int spPhIndex = 0; spPhIndex < brdf->getNumSpecPhi();   ++spPhIndex) {
                const Spectrum& minSp = brdf->getSpectrum(inThIndex, minInPhiIndex, spThIndex, spPhIndex);
                const Spectrum& maxSp = brdf->getSpectrum(inThIndex, maxInPhiIndex, spThIndex, spPhIndex);
                Spectrum sp = (minSp + maxSp) / 2.0;

                brdf->setSpectrum(inThIndex, minInPhiIndex, spThIndex, spPhIndex, sp);
                brdf->setSpectrum(inThIndex, maxInPhiIndex, spThIndex, spPhIndex, sp);
            }}}
        }
    }

    if (brdf->getNumSpecPhi() >= 2) {
        if (isEqual(brdf->getSpecTheta(0), 0.0f)) {
            // Equalize samples for specular azimuthal angles if a specular polar angle is 0.
            for (int inThIndex = 0; inThIndex < brdf->getNumInTheta();   ++inThIndex) {
            for (int inPhIndex = 0; inPhIndex < brdf->getNumInPhi();     ++inPhIndex) {
                Arrayd sumSp = Arrayd::Zero(ss->getNumWavelengths());
                for (int spPhIndex = 0; spPhIndex < brdf->getNumSpecPhi(); ++spPhIndex) {
                    sumSp += brdf->getSpectrum(inThIndex, inPhIndex, 0, spPhIndex).cast<Arrayd::Scalar>();
                }

                Spectrum sp = sumSp.cast<Spectrum::Scalar>() / brdf->getNumSpecPhi();
                for (int spPhIndex = 0; spPhIndex < brdf->getNumSpecPhi(); ++spPhIndex) {
                    brdf->setSpectrum(inThIndex, inPhIndex, 0, spPhIndex, sp);
                }
            }}
        }

        if (isEqual(brdf->getSpecTheta(brdf->getNumSpecTheta() - 1), SpecularCoordinateSystem::MAX_ANGLE2)) {
            // Equalize samples for specular azimuthal angles if a specular polar angle is 2PI.
            for (int inThIndex = 0; inThIndex < brdf->getNumInTheta(); ++inThIndex) {
            for (int inPhIndex = 0; inPhIndex < brdf->getNumInPhi();   ++inPhIndex) {
                Arrayd sumSp = Arrayd::Zero(ss->getNumWavelengths());
                for (int spPhIndex = 0; spPhIndex < brdf->getNumSpecPhi(); ++spPhIndex) {
                    sumSp += brdf->getSpectrum(inThIndex, inPhIndex, brdf->getNumSpecTheta() - 1, spPhIndex).cast<Arrayd::Scalar>();
                }

                Spectrum sp = sumSp.cast<Spectrum::Scalar>() / brdf->getNumSpecPhi();
                for (int spPhIndex = 0; spPhIndex < brdf->getNumSpecPhi(); ++spPhIndex) {
                    brdf->setSpectrum(inThIndex, inPhIndex, brdf->getNumSpecTheta() - 1, spPhIndex, sp);
                }
            }}
        }

        int minSpPhiIndex = 0;
        int maxSpPhiIndex = brdf->getNumSpecPhi() - 1;
        if (isEqual(brdf->getSpecPhi(minSpPhiIndex), SpecularCoordinateSystem::MIN_ANGLE3) &&
            isEqual(brdf->getSpecPhi(maxSpPhiIndex), SpecularCoordinateSystem::MAX_ANGLE3)) {
            // Equalize samples with overlapping specular azimuthal angles.
            for (int inThIndex = 0; inThIndex < brdf->getNumInTheta();   ++inThIndex) {
            for (int inPhIndex = 0; inPhIndex < brdf->getNumInPhi();     ++inPhIndex) {
            for (int spThIndex = 0; spThIndex < brdf->getNumSpecTheta(); ++spThIndex) {
                const Spectrum& minSp = brdf->getSpectrum(inThIndex, inPhIndex, spThIndex, minSpPhiIndex);
                const Spectrum& maxSp = brdf->getSpectrum(inThIndex, inPhIndex, spThIndex, maxSpPhiIndex);
                Spectrum sp = (minSp + maxSp) / 2.0;

                brdf->setSpectrum(inThIndex, inPhIndex, spThIndex, minSpPhiIndex, sp);
                brdf->setSpectrum(inThIndex, inPhIndex, spThIndex, maxSpPhiIndex, sp);
            }}}
        }
    }
}

void lb::removeSpecularValues(SpecularCoordinatesBrdf* brdf, float maxSpecularTheta)
{
    if (brdf->getNumSpecTheta() == 1 ||
        brdf->getNumSpecPhi() == 1) {
        return;
    }

    const SampleSet* ss = brdf->getSampleSet();

    int spThBoundary = 0;

    // Extrapolate values in specular directions and search the index of boundaries.
    for (int inThIndex = 0; inThIndex < brdf->getNumInTheta();       ++inThIndex) {
    for (int inPhIndex = 0; inPhIndex < brdf->getNumInPhi();         ++inPhIndex) {
    for (int spThIndex = 0; spThIndex < brdf->getNumSpecTheta() - 1; ++spThIndex) {
        if (brdf->getSpecTheta(spThIndex) <= maxSpecularTheta) {
            continue;
        }

        Spectrum sumSp = Spectrum::Zero(ss->getNumWavelengths());

        int numSpectra = 0;

        for (int spPhIndex = 0; spPhIndex < brdf->getNumSpecPhi() - 1; ++spPhIndex) {
            Vec3 inDir, outDir;
            brdf->toXyz(brdf->getInTheta(inThIndex),
                        brdf->getInPhi(inPhIndex),
                        brdf->getSpecTheta(spThIndex),
                        brdf->getSpecPhi(spPhIndex),
                        &inDir, &outDir);

            if (isDownwardDir(outDir)) {
                continue;
            }

            Spectrum sp     = brdf->getSpectrum(inThIndex, inPhIndex, spThIndex, spPhIndex);
            Spectrum nextSp = brdf->getSpectrum(inThIndex, inPhIndex, spThIndex + 1, spPhIndex);

            float spTh     = brdf->getSpecTheta(spThIndex);
            float nextSpTh = brdf->getSpecTheta(spThIndex + 1);

            float ratio = (0.0f - spTh) / (nextSpTh - spTh);
            Spectrum extrapolatedSp = lerp(sp, nextSp, ratio);
            sumSp += extrapolatedSp.cwiseMax(0.0);

            ++numSpectra;
        }

        if (numSpectra == 0) {
            break;
        }

        Spectrum specularSp = sumSp / numSpectra;

        for (int spPhIndex = 0; spPhIndex < brdf->getNumSpecPhi(); ++spPhIndex) {
            brdf->setSpectrum(inThIndex, inPhIndex, 0, spPhIndex, specularSp);
        }

        spThBoundary = spThIndex;

        break;
    }}}

    // Interpolate values between specular directions and the maximum specular polar angle.
    for (int inThIndex = 0; inThIndex < brdf->getNumInTheta();   ++inThIndex) {
    for (int inPhIndex = 0; inPhIndex < brdf->getNumInPhi();     ++inPhIndex) {
    for (int spThIndex = 1; spThIndex < brdf->getNumSpecTheta(); ++spThIndex) {
        if (brdf->getSpecTheta(spThIndex) > maxSpecularTheta) {
            break;
        }

        Spectrum specularSp = brdf->getSpectrum(inThIndex, inPhIndex, 0, 0);

        for (int spPhIndex = 0; spPhIndex < brdf->getNumSpecPhi(); ++spPhIndex) {
            Spectrum boundarySp = brdf->getSpectrum(inThIndex, inPhIndex, spThBoundary, spPhIndex);
            float boundarySpTh = brdf->getSpecTheta(spThBoundary);

            float ratio = brdf->getSpecTheta(spThIndex) / boundarySpTh;

            // Smooth a peak.
            if (ratio < 0.5f) {
                ratio = smoothstep(0.0f, boundarySpTh, brdf->getSpecTheta(spThIndex));
            }

            Spectrum sp = lerp(specularSp, boundarySp, ratio);
            brdf->setSpectrum(inThIndex, inPhIndex, spThIndex, spPhIndex, sp);
        }
    }}}
}

/* Insert a BRDF in a base BRDF along incoming azimuthal angle. */
template <typename BrdfT>
BrdfT* insertBrdfAlongInPhiTemplate(const BrdfT&    baseBrdf,
                                    const BrdfT&    insertedBrdf,
                                    float           inPhi)
{
    const SampleSet* baseSs     = baseBrdf.getSampleSet();
    const SampleSet* insertedSs = insertedBrdf.getSampleSet();

    if (!hasSameColor(*baseSs, *insertedSs)) {
        lbError << "[lb::insertBrdfAlongInPhiTemplate] Color models or wavelengths do not match.";
        return 0;
    }

    if (insertedSs->getNumAngles1() != 1) {
        lbError
            << "[lb::insertBrdfAlongInPhi] The number of incoming azimuthal angles must be 1. The number of angles: "
            << insertedSs->getNumAngles1();
        return 0;
    }

    if (inPhi < 0.0f || inPhi > 2.0f * PI_F) {
        lbError
            << "[lb::insertBrdfAlongInPhi] Specified incoming azimuthal angle is out of range: "
            << inPhi;
        return 0;
    }

    // Find the index of the inserted angle.
    int insertedIndex = baseSs->getNumAngles1();
    for (int i = 0; i < baseSs->getNumAngles1(); ++i) {
        if (isEqual(baseSs->getAngle1(i), inPhi)) {
            lbError
                << "[lb::insertBrdfAlongInPhi] Specified incoming azimuthal angle is already used: "
                << inPhi;
            return 0;
        }

        if (baseSs->getAngle1(i) > inPhi) {
            insertedIndex = i;
            break;
        }
    }

    BrdfT* brdf = baseBrdf.clone();
    SampleSet* ss = brdf->getSampleSet();

    ss->resizeAngles(baseSs->getNumAngles0(),
                     baseSs->getNumAngles1() + 1,
                     baseSs->getNumAngles2(),
                     baseSs->getNumAngles3());

    ss->getAngles0() = baseSs->getAngles0();
    ss->getAngles2() = baseSs->getAngles2();
    ss->getAngles3() = baseSs->getAngles3();

    // Add the inserted angle.
    Arrayf& angles1 = ss->getAngles1();
    for (int i = 0; i < baseSs->getNumAngles1(); ++i) {
        angles1[i] = baseSs->getAngles1()[i];
    }
    angles1[angles1.size() - 1] = inPhi;
    std::sort(angles1.data(), angles1.data() + angles1.size());

    ss->updateAngleAttributes();

    for (int i0 = 0; i0 < ss->getNumAngles0(); ++i0) {
    for (int i1 = 0; i1 < ss->getNumAngles1(); ++i1) {
    for (int i2 = 0; i2 < ss->getNumAngles2(); ++i2) {
    for (int i3 = 0; i3 < ss->getNumAngles3(); ++i3) {
        Spectrum sp;
        if (i1 < insertedIndex) {
            sp = baseSs->getSpectrum(i0, i1, i2, i3);
        }
        else if (i1 == insertedIndex) {
            float inTheta = ss->getAngle0(i0);

            // An incoming polar angle of zero is offset to validate an incoming azimuthal angle.
            float offsetInTheta = std::max(inTheta, EPSILON_F);

            Vec3 inDir, outDir;
            brdf->toXyz(offsetInTheta,
                        ss->getAngle1(i1),
                        ss->getAngle2(i2),
                        ss->getAngle3(i3),
                        &inDir, &outDir);
            sp = insertedBrdf.getSpectrum(inDir, outDir);
        }
        else {
            sp = baseSs->getSpectrum(i0, i1 - 1, i2, i3);
        }

        ss->setSpectrum(i0, i1, i2, i3, sp);
    }}}}

    return brdf;
}

SpecularCoordinatesBrdf* lb::insertBrdfAlongInPhi(const SpecularCoordinatesBrdf&    baseBrdf,
                                                  const SpecularCoordinatesBrdf&    insertedBrdf,
                                                  float                             inPhi)
{
    return insertBrdfAlongInPhiTemplate(baseBrdf, insertedBrdf, inPhi);
}

SphericalCoordinatesBrdf* lb::insertBrdfAlongInPhi(const SphericalCoordinatesBrdf&  baseBrdf,
                                                   const SphericalCoordinatesBrdf&  insertedBrdf,
                                                   float                            inPhi)
{
    return insertBrdfAlongInPhiTemplate(baseBrdf, insertedBrdf, inPhi);
}

void lb::extrapolateSamplesWithReflectances(SpecularCoordinatesBrdf* brdf, float incomingTheta)
{
    if (brdf->getNumInTheta() < 3 ||
        brdf->getInTheta(1) > incomingTheta ||
        brdf->getInTheta(brdf->getNumInTheta() - 1) <= incomingTheta) {
        return;
    }

    Spectrum diffuseThresholds = findDiffuseThresholds(*brdf, incomingTheta);

    SpecularCoordinatesBrdf* glossyBrdf = brdf->clone();
    SpecularCoordinatesBrdf* diffuseBrdf = brdf->clone();

    editComponents(*brdf, glossyBrdf,  diffuseThresholds, 1.0f, 1.0f, 0.0f);
    editComponents(*brdf, diffuseBrdf, diffuseThresholds, 0.0f, 1.0f, 1.0f);

    SampleSet2D* gRefs = computeReflectances(*glossyBrdf);
    SampleSet2D* dRefs = computeReflectances(*diffuseBrdf);

    int inThBoundaryIndex = 0;
    for (int inThIndex = 0; inThIndex < brdf->getNumInTheta(); ++inThIndex) {
        if (brdf->getInTheta(inThIndex) > incomingTheta) {
            break;
        }

        inThBoundaryIndex = inThIndex;
    }

    for (int inThIndex = inThBoundaryIndex + 1; inThIndex < brdf->getNumInTheta();   ++inThIndex) {
    for (int inPhIndex = 0;                     inPhIndex < brdf->getNumInPhi();     ++inPhIndex) {
    for (int spThIndex = 0;                     spThIndex < brdf->getNumSpecTheta(); ++spThIndex) {
        Spectrum gRef0 = gRefs->getSpectrum(inThBoundaryIndex - 1, inPhIndex);
        Spectrum gRef1 = gRefs->getSpectrum(inThBoundaryIndex    , inPhIndex);

        Spectrum dRef0 = dRefs->getSpectrum(inThBoundaryIndex - 1, inPhIndex);
        Spectrum dRef1 = dRefs->getSpectrum(inThBoundaryIndex,     inPhIndex);

        float angle0 = brdf->getInTheta(inThBoundaryIndex - 1);
        float angle1 = brdf->getInTheta(inThBoundaryIndex);
        float t = (brdf->getInTheta(inThIndex) - angle0) / (angle1 - angle0);

        Spectrum extrapolatedGRef = lerp(gRef0, gRef1, t);
        Spectrum extrapolatedDRef = lerp(dRef0, dRef1, t);

        Spectrum gRef = gRefs->getSpectrum(inThIndex, inPhIndex);
        Spectrum dRef = dRefs->getSpectrum(inThIndex, inPhIndex);

        // Avoid dividing by zero.
        const float minRef = 0.01f;
        extrapolatedGRef = extrapolatedGRef.cwiseMax(minRef);
        extrapolatedDRef = extrapolatedDRef.cwiseMax(minRef);
        gRef = gRef.cwiseMax(minRef);
        dRef = dRef.cwiseMax(minRef);

        for (int spPhIndex = 0; spPhIndex < brdf->getNumSpecPhi(); ++spPhIndex) {
            Spectrum gSp = glossyBrdf->getSpectrum(inThIndex, inPhIndex, spThIndex, spPhIndex);
            Spectrum dSp = diffuseBrdf->getSpectrum(inThIndex, inPhIndex, spThIndex, spPhIndex);

            gSp *= extrapolatedGRef / gRef;
            dSp *= extrapolatedDRef / dRef;

            Spectrum gdSp = (gSp + dSp).cwiseMax(0.0f);
            brdf->setSpectrum(inThIndex, inPhIndex, spThIndex, spPhIndex, gdSp);
        }
    }}}

    delete diffuseBrdf;
    delete glossyBrdf;
    delete dRefs;
    delete gRefs;
}

void lb::copySpectraFromPhiOf0To360(SampleSet* samples)
{
    if (samples->getNumAngles1() >= 2 &&
        samples->getAngle1(0) == 0.0f &&
        samples->getAngle1(samples->getNumAngles1() - 1) >= SphericalCoordinateSystem::MAX_ANGLE1) {
        for (int i0 = 0; i0 < samples->getNumAngles0(); ++i0) {
        for (int i2 = 0; i2 < samples->getNumAngles2(); ++i2) {
        for (int i3 = 0; i3 < samples->getNumAngles3(); ++i3) {
            const Spectrum& sp = samples->getSpectrum(i0, 0, i2, i3);
            samples->setSpectrum(i0, samples->getNumAngles1() - 1, i2, i3, sp);
        }}}
    }

    if (samples->getNumAngles3() >= 2 &&
        samples->getAngle3(0) == 0.0f &&
        samples->getAngle3(samples->getNumAngles3() - 1) >= SphericalCoordinateSystem::MAX_ANGLE3) {
        for (int i0 = 0; i0 < samples->getNumAngles0(); ++i0) {
        for (int i1 = 0; i1 < samples->getNumAngles1(); ++i1) {
        for (int i2 = 0; i2 < samples->getNumAngles2(); ++i2) {
            const Spectrum& sp = samples->getSpectrum(i0, i1, i2, 0);
            samples->setSpectrum(i0, i1, i2, samples->getNumAngles3() - 1, sp);
        }}}
    }
}

bool lb::fillSpectraAtInThetaOf90(Brdf* brdf, Spectrum::Scalar value)
{
    if (!dynamic_cast<SpecularCoordinatesBrdf*>(brdf) &&
        !dynamic_cast<SphericalCoordinatesBrdf*>(brdf)) {
        lbError << "[fillSpectraAtInThetaOf90] Unsupported type of BRDF";
        return false;
    }

    SampleSet* ss = brdf->getSampleSet();

    int endIndex0 = ss->getNumAngles0() - 1;
    if (!isEqual(ss->getAngle0(endIndex0), PI_2_F)) {
        return false;
    }

    for (int i1 = 0; i1 < ss->getNumAngles1(); ++i1) {
    for (int i2 = 0; i2 < ss->getNumAngles2(); ++i2) {
    for (int i3 = 0; i3 < ss->getNumAngles3(); ++i3) {
        Spectrum& sp = ss->getSpectrum(endIndex0, i1, i2, i3);
        sp.fill(value);
    }}}

    return true;
}

void lb::xyzToSrgb(SampleSet* samples)
{
    ColorModel cm = samples->getColorModel();
    if (cm != XYZ_MODEL) {
        lbError << "[xyzToSrgb] Not CIE-XYZ model: " << cm;
        return;
    }

    for (int i0 = 0; i0 < samples->getNumAngles0(); ++i0) {
    for (int i1 = 0; i1 < samples->getNumAngles1(); ++i1) {
    for (int i2 = 0; i2 < samples->getNumAngles2(); ++i2) {
    for (int i3 = 0; i3 < samples->getNumAngles3(); ++i3) {
        const Spectrum& xyz = samples->getSpectrum(i0, i1, i2, i3);
        Spectrum rgb = xyzToSrgb<Vec3f>(xyz);
        samples->setSpectrum(i0, i1, i2, i3, rgb);
    }}}}

    samples->setColorModel(RGB_MODEL);
}

void lb::fillSpectra(SampleSet* samples, Spectrum::Scalar value)
{
    fillSpectra(samples->getSpectra(), value);
}

void lb::fillSpectra(SpectrumList& spectra, Spectrum::Scalar value)
{
    for (auto& sp : spectra) {
        sp.fill(value);
    }
}

bool lb::subtract(const Brdf& src0, const Brdf& src1, Brdf* dest)
{
    const SampleSet* ss0 = src0.getSampleSet();
    const SampleSet* ss1 = src1.getSampleSet();
    SampleSet* ss = dest->getSampleSet();

    if (ss0->getColorModel() != ss1->getColorModel() ||
        ss0->getColorModel() != ss->getColorModel()) {
        lbError << "[subtract] Color models are not identical.";
        return false;
    }

    const Arrayf& wls0 = ss0->getWavelengths();
    const Arrayf& wls1 = ss1->getWavelengths();
    const Arrayf& wls = ss->getWavelengths();
    if (!wls0.isApprox(wls1) ||
        !wls0.isApprox(wls)) {
        lbError << "[subtract] Wavelengths are not identical.";
        return false;
    }

    for (int i0 = 0; i0 < ss->getNumAngles0(); ++i0) {
    for (int i1 = 0; i1 < ss->getNumAngles1(); ++i1) {
    for (int i2 = 0; i2 < ss->getNumAngles2(); ++i2) {
    for (int i3 = 0; i3 < ss->getNumAngles3(); ++i3) {
        Vec3 inDir, outDir;
        dest->getInOutDirection(i0, i1, i2, i3, &inDir, &outDir);

        const Spectrum& sp0 = src0.getSpectrum(inDir, outDir);
        const Spectrum& sp1 = src1.getSpectrum(inDir, outDir);

        ss->setSpectrum(i0, i1, i2, i3, sp0 - sp1);
    }}}}

    return true;
}

void lb::multiplySpectra(SampleSet* samples, Spectrum::Scalar value)
{
    for (auto& sp : samples->getSpectra()) {
        sp *= value;
    }
}

void lb::fixNegativeSpectra(SampleSet* samples)
{
    fixNegativeSpectra(samples->getSpectra());
}

void lb::fixNegativeSpectra(SpectrumList& spectra)
{
    for (auto& sp : spectra) {
        sp = sp.cwiseMax(0);
    }
}
