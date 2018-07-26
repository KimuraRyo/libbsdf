// =================================================================== //
// Copyright (C) 2014-2018 Kimura Ryo                                  //
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

#include <libbsdf/Brdf/Brdf.h>
#include <libbsdf/Brdf/SpecularCoordinatesBrdf.h>
#include <libbsdf/Brdf/SphericalCoordinatesBrdf.h>

#include <libbsdf/Common/PoissonDiskDistributionOnSphere.h>
#include <libbsdf/Common/SpectrumUtility.h>
#include <libbsdf/Common/SphericalCoordinateSystem.h>

#include <libbsdf/ReflectanceModel/Fresnel.h>

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

        Spectrum& sp = ss->getSpectrum(i0, i1, i2, i3);

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

    Integrator integrator(PoissonDiskDistributionOnSphere::NUM_SAMPLES_ON_HEMISPHERE, true);

    for (int inThIndex = 0; inThIndex < brdf->getNumInTheta(); ++inThIndex) {
    for (int inPhIndex = 0; inPhIndex < brdf->getNumInPhi();   ++inPhIndex) {
        Vec3 inDir = SphericalCoordinateSystem::toXyz(brdf->getInTheta(inThIndex),
                                                      brdf->getInPhi(inPhIndex));
        Spectrum sp = integrator.computeReflectance(*brdf, inDir);

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

void lb::fixEnergyConservation(SpecularCoordinatesBrdf* brdf,
                               const SampleSet2D&       specularReflectances)
{
    SampleSet* ss = brdf->getSampleSet();

    Integrator integrator(PoissonDiskDistributionOnSphere::NUM_SAMPLES_ON_HEMISPHERE, true);

    for (int inThIndex = 0; inThIndex < brdf->getNumInTheta(); ++inThIndex) {
    for (int inPhIndex = 0; inPhIndex < brdf->getNumInPhi();   ++inPhIndex) {
        Vec3 inDir = SphericalCoordinateSystem::toXyz(brdf->getInTheta(inThIndex),
                                                      brdf->getInPhi(inPhIndex));
        Spectrum sp = integrator.computeReflectance(*brdf, inDir);

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

        Spectrum specRefSp = specularReflectances.getSpectrum(inDir);
        sp += specRefSp;

        maxReflectance = 0;
        int maxIndex;
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
                const float coeff = 0.999546f; // Reflectance of Lambertian using lb::Integrator.
                fixedSp /= maxReflectance / coeff;
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
                const float coeff = 0.999546f; // Reflectance of Lambertian using lb::Integrator.
                fixedSp /= maxReflectance / coeff;
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
        bool found = false;
        
        int boundary0;
        // Search a boundary.
        for (int spPhIndex = 0; spPhIndex < brdf->getNumSpecPhi(); ++spPhIndex) {
            Vec3 inDir, outDir;
            brdf->toXyz(brdf->getInTheta(inThIndex),
                        brdf->getInPhi(inPhIndex),
                        brdf->getSpecTheta(spThIndex),
                        brdf->getSpecPhi(spPhIndex),
                        &inDir, &outDir);

            boundary0 = spPhIndex;
            if (outDir.z() < 0.0) {
                found = true;
                break;
            }
        }

        if (!found) {
            continue;
        }

        int boundary1;
        // Search a boundary from the opposite direction.
        for (int spPhIndex = brdf->getNumSpecPhi() - 1; spPhIndex >= 0; --spPhIndex) {
            Vec3 inDir, outDir;
            brdf->toXyz(brdf->getInTheta(inThIndex),
                        brdf->getInPhi(inPhIndex),
                        brdf->getSpecTheta(spThIndex),
                        brdf->getSpecPhi(spPhIndex),
                        &inDir, &outDir);

            boundary1 = spPhIndex;
            if (outDir.z() < 0.0) {
                break;
            }
        }

        // Fill values.
        for (int spPhIndex = 0; spPhIndex < brdf->getNumSpecPhi(); ++spPhIndex) {
            Vec3 inDir, outDir;
            brdf->toXyz(brdf->getInTheta(inThIndex),
                        brdf->getInPhi(inPhIndex),
                        brdf->getSpecTheta(spThIndex),
                        brdf->getSpecPhi(spPhIndex),
                        &inDir, &outDir);

            if (outDir.z() >= 0.0) {
                continue;
            }

            float spPh = brdf->getSpecPhi(spPhIndex);
            float spPh0 = brdf->getSpecPhi(boundary0);
            float spPh1 = brdf->getSpecPhi(boundary1);

            int boundary;
            if (spPh - spPh0 <= spPh1 - spPh) {
                boundary = boundary0;
            }
            else {
                boundary = boundary1;
            }

            Spectrum sp = brdf->getSpectrum(inThIndex, inPhIndex, spThIndex, boundary);
            brdf->setSpectrum(inThIndex, inPhIndex, spThIndex, spPhIndex, sp);
        }
    }}}
}

void lb::removeSpecularValues(SpecularCoordinatesBrdf* brdf, float maxSpecularTheta)
{
    if (brdf->getNumSpecTheta() == 1 ||
        brdf->getNumSpecPhi() == 1) {
        return;
    }

    const SampleSet* ss = brdf->getSampleSet();

    int spThBoundary;

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

            if (outDir.z() < 0.0) {
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

Brdf* lb::insertBrdfAlongInPhi(const SphericalCoordinatesBrdf&  baseBrdf,
                               const SphericalCoordinatesBrdf&  insertedBrdf,
                               float                            inPhi)
{
    //bool inPhiUsed = (dynamic_cast<const SphericalCoordinatesBrdf*>(&baseBrdf) &&
    //                  dynamic_cast<const SphericalCoordinatesBrdf*>(&insertedBrdf)) ||
    //                 (dynamic_cast<const SpecularCoordinatesBrdf*>(&baseBrdf) &&
    //                  dynamic_cast<const SpecularCoordinatesBrdf*>(&insertedBrdf));

    //if (!inPhiUsed) {
    //    std::cerr
    //        << "[lb::insertBrdfAlongInPhi] Invalid coordinate system is used."
    //        << std::endl;
    //    return 0;
    //}

    const SampleSet* baseSs = baseBrdf.getSampleSet();
    const SampleSet* insertedSs = insertedBrdf.getSampleSet();

    //if (!hasSameColor(baseSs, insertedSs)) {
    //    std::cerr
    //        << "[lb::insertBrdfAlongInPhi] Color models or wavelengths do not match."
    //        << std::endl;
    //    return 0;
    //}

    if (baseSs->getColorModel() != insertedSs->getColorModel()) {
        std::cerr
            << "[Brdf::insertBrdfAlongInPhi] Color models do not match: "
            << baseSs->getColorModel() << ", " << insertedSs->getColorModel()
            << std::endl;
        return 0;
    }

    if (baseSs->getNumWavelengths() != insertedSs->getNumWavelengths() ||
        !baseSs->getWavelengths().isApprox(insertedSs->getWavelengths())) {
        std::cerr
            << "[lb::insertBrdfAlongInPhi] Wavelengths do not match: "
            << baseSs->getWavelengths() << ", " << insertedSs->getWavelengths()
            << std::endl;
        return 0;
    }

    if (insertedSs->getNumAngles1() != 1) {
        std::cerr
            << "[lb::insertBrdfAlongInPhi] The number of incoming azimuthal angles must be 1. insertedSs->getNumAngles1(): "
            << insertedSs->getNumAngles1() << std::endl;
        return 0;
    }

    if (inPhi < 0.0f || inPhi > 2.0 * PI_F) {
        std::cerr << "[lb::insertBrdfAlongInPhi] inPhi is out of range: " << inPhi << std::endl;
        return 0;
    }

    // Find the index of the inserted angle.
    int insertedIndex = baseSs->getNumAngles1();
    for (int i = 0; i < baseSs->getNumAngles1(); ++i) {
        if (baseSs->getAngle1(i) == inPhi) {
            std::cerr << "[lb::insertBrdfAlongInPhi] inPhi is already used: " << inPhi << std::endl;
            return 0;
        }

        if (baseSs->getAngle1(i) > inPhi) {
            insertedIndex = i;
            break;
        }
    }

    Brdf* brdf = baseBrdf.clone();
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

    for (int i0 = 0; i0 < ss->getNumAngles0(); ++i0) {
    for (int i1 = 0; i1 < ss->getNumAngles1(); ++i1) {
    for (int i2 = 0; i2 < ss->getNumAngles2(); ++i2) {
    for (int i3 = 0; i3 < ss->getNumAngles3(); ++i3) {
        Spectrum sp;
        if (i1 < insertedIndex) {
            sp = baseSs->getSpectrum(i0, i1, i2, i3);
        }
        else if (i1 == insertedIndex) {
            Vec3 inDir, outDir;
            insertedBrdf.toXyz(insertedSs->getAngle0(i0),
                               insertedSs->getAngle1(i1),
                               insertedSs->getAngle2(i2),
                               insertedSs->getAngle3(i3) - inPhi,
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

void lb::copySpectraFromPhiOfZeroTo2PI(SampleSet* samples)
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
        std::cerr << "[fillSpectraAtInThetaOf90] Unsupported type of BRDF" << std::endl;
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
        std::cerr << "[xyzToSrgb] Not CIE-XYZ model: " << cm << std::endl;
        return;
    }

    for (int i0 = 0; i0 < samples->getNumAngles0(); ++i0) {
    for (int i1 = 0; i1 < samples->getNumAngles1(); ++i1) {
    for (int i2 = 0; i2 < samples->getNumAngles2(); ++i2) {
    for (int i3 = 0; i3 < samples->getNumAngles3(); ++i3) {
        const Spectrum& xyz = samples->getSpectrum(i0, i1, i2, i3);
        Spectrum rgb = xyzToSrgb(xyz);
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
    for (auto it = spectra.begin(); it != spectra.end(); ++it) {
        it->fill(value);
    }
}

bool lb::subtract(const Brdf& src0, const Brdf& src1, Brdf* dest)
{
    const SampleSet* ss0 = src0.getSampleSet();
    const SampleSet* ss1 = src1.getSampleSet();
    SampleSet* ss = dest->getSampleSet();

    if (ss0->getColorModel() != ss1->getColorModel() ||
        ss0->getColorModel() != ss->getColorModel()) {
        std::cerr << "[subtract] Color models are not identical." << std::endl;
        return false;
    }

    const Arrayf& wls0 = ss0->getWavelengths();
    const Arrayf& wls1 = ss1->getWavelengths();
    const Arrayf& wls = ss->getWavelengths();
    if (!wls0.isApprox(wls1) ||
        !wls0.isApprox(wls)) {
        std::cerr << "[subtract] Wavelengths are not identical." << std::endl;
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
    SpectrumList& spectra = samples->getSpectra();
    for (auto it = spectra.begin(); it != spectra.end(); ++it) {
        *it *= value;
    }
}

void lb::fixNegativeSpectra(SampleSet* samples)
{
    fixNegativeSpectra(samples->getSpectra());
}

void lb::fixNegativeSpectra(SpectrumList& spectra)
{
    for (auto it = spectra.begin(); it != spectra.end(); ++it) {
        *it = it->cwiseMax(0.0f);
    }
}
