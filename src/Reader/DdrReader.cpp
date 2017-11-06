// =================================================================== //
// Copyright (C) 2014-2017 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#include <libbsdf/Reader/DdrReader.h>

#include <fstream>

#include <libbsdf/Brdf/Integrator.h>
#include <libbsdf/Common/PoissonDiskDistributionOnSphere.h>
#include <libbsdf/Reader/DdrSdrUtility.h>

using namespace lb;

SpecularCoordinatesBrdf* DdrReader::read(const std::string& fileName)
{
    // std::ios_base::binary is used to read line endings of CR+LF and LF.
    std::ifstream ifs(fileName.c_str(), std::ios_base::binary);
    if (ifs.fail()) {
        std::cerr << "[DdrReader::read] Could not open: " << fileName << std::endl;
        return 0;
    }

    std::ios_base::sync_with_stdio(false);

    SourceType sourceType = UNKNOWN_SOURCE;

    ddr_sdr_utility::SymmetryType symmetryType  = ddr_sdr_utility::PLANE_SYMMETRICAL;
    ddr_sdr_utility::UnitType     unitType      = ddr_sdr_utility::LUMINANCE_ABSOLUTE;

    ColorModel colorModel = RGB_MODEL;

    int numWavelengths = 1;

    std::vector<float> inThetaDegrees;
    std::vector<float> inPhiDegrees;
    std::vector<float> spThetaDegrees;
    std::vector<float> spPhiDegrees;

    ddr_sdr_utility::ignoreCommentLines(ifs);
    std::ifstream::pos_type pos = ifs.tellg();

    using reader_utility::isEqual;

    // Read a header.
    std::string headStr;
    while (ifs >> headStr) {
        ddr_sdr_utility::ignoreCommentLines(ifs);

        if (headStr.empty()) {
            continue;
        }
        else if (isEqual(headStr, "Source")) {
            std::string typeStr;
            ifs >> typeStr;

            if (isEqual(typeStr, "Measured")) {
                sourceType = MEASURED_SOURCE;
            }
            else if (isEqual(typeStr, "Generated")) {
                sourceType = GENERATED_SOURCE;
            }
            else if (isEqual(typeStr, "Edited")) {
                sourceType = EDITED_SOURCE;
            }
            else if (isEqual(typeStr, "Morphed")) {
                sourceType = UNKNOWN_SOURCE;
            }
            else {
                reader_utility::logNotImplementedKeyword(typeStr);
            }
        }
        else if (isEqual(headStr, "TypeSym")) {
            std::string typeStr;
            ifs >> typeStr;

            if (isEqual(typeStr, "AxiSymmetrical")) {
                symmetryType = ddr_sdr_utility::AXI_SYMMETRICAL;
            }
            else if (isEqual(typeStr, "DirSymmetrical")) {
                symmetryType = ddr_sdr_utility::DIR_SYMMETRICAL;
            }
            else if (isEqual(typeStr, "PlaneSymmetrical")) {
                symmetryType = ddr_sdr_utility::PLANE_SYMMETRICAL;
            }
            else if (isEqual(typeStr, "ASymmetrical")) {
                std::ifstream::pos_type asymmetricalPos = ifs.tellg();
                ifs >> typeStr;

                if (isEqual(typeStr, "4D")) {
                    symmetryType = ddr_sdr_utility::ASYMMETRICAL_4D;
                }
                else {
                    symmetryType = ddr_sdr_utility::ASYMMETRICAL;
                    ifs.seekg(asymmetricalPos, std::ios_base::beg);
                }
            }
            else {
                reader_utility::logNotImplementedKeyword(typeStr);
                return 0;
            }
        }
        else if (isEqual(headStr, "TypeColorModel")) {
            std::string typeStr;
            ifs >> typeStr;

            if (isEqual(typeStr, "rgb")) {
                colorModel = RGB_MODEL;
                numWavelengths = 3;
            }
            else if (isEqual(typeStr, "spectral")) {
                colorModel = SPECTRAL_MODEL;
                ifs >> numWavelengths;
            }
            else if (isEqual(typeStr, "bw")) {
                colorModel = MONOCHROMATIC_MODEL;
                numWavelengths = 1;
            }
            else {
                reader_utility::logNotImplementedKeyword(typeStr);
                return 0;
            }
        }
        else if (isEqual(headStr, "TypeData")) {
            std::string typeStr;
            ifs >> typeStr;

            if (isEqual(typeStr, "Luminance")) {
                ifs >> typeStr;

                if (isEqual(typeStr, "Absolute")) {
                    unitType = ddr_sdr_utility::LUMINANCE_ABSOLUTE;
                }
                else if (isEqual(typeStr, "Relative")) {
                    unitType = ddr_sdr_utility::LUMINANCE_RELATIVE;
                }
                else {
                    reader_utility::ignoreLine(ifs);
                }
            }
            else if (isEqual(typeStr, "Intensity")) {
                ifs >> typeStr;

                if (isEqual(typeStr, "Absolute")) {
                    unitType = ddr_sdr_utility::INTENSITY_ABSOLUTE;
                }
                else if (isEqual(typeStr, "Relative")) {
                    unitType = ddr_sdr_utility::INTENSITY_RELATIVE;
                }
                else {
                    reader_utility::ignoreLine(ifs);
                }
                reader_utility::ignoreLine(ifs);
            }
            else {
                reader_utility::logNotImplementedKeyword(typeStr);
                return 0;
            }
        }
        else if (isEqual(headStr, "psi")) {
            int numInPhi;
            ifs >> numInPhi;
            for (int i = 0; i < numInPhi; ++i) {
                float angle;
                ifs >> angle;
                inPhiDegrees.push_back(angle);
            }
        }
        else if (isEqual(headStr, "sigma")) {
            int numInTheta;
            ifs >> numInTheta;
            for (int i = 0; i < numInTheta; ++i) {
                float angle;
                ifs >> angle;
                inThetaDegrees.push_back(angle);
            }
        }
        else if (isEqual(headStr, "sigmaT")) {
            reader_utility::logNotImplementedKeyword(headStr);
            reader_utility::ignoreLine(ifs);
        }
        else if (isEqual(headStr, "phi")) {
            int numSpecPhi;
            ifs >> numSpecPhi;
            for (int i = 0; i < numSpecPhi; ++i) {
                float angle;
                ifs >> angle;
                spPhiDegrees.push_back(angle);
            }
        }
        else if (isEqual(headStr, "theta")) {
            int numSpecTheta;
            ifs >> numSpecTheta;
            for (int i = 0; i < numSpecTheta; ++i) {
                float angle;
                ifs >> angle;
                spThetaDegrees.push_back(angle);
            }
        }
        else if (isEqual(headStr, "wl") ||
                 isEqual(headStr, "bw") ||
                 isEqual(headStr, "red") ||
                 isEqual(headStr, "gre") ||
                 isEqual(headStr, "blu") ||
                 isEqual(headStr, "green") ||   // Not correct specification
                 isEqual(headStr, "blue")) {    // Not correct specification
            ifs.seekg(pos, std::ios_base::beg);
            break;
        }

        pos = ifs.tellg();
    }

    if (inThetaDegrees.empty() ||
        spThetaDegrees.empty() ||
        spPhiDegrees.empty()) {
        std::cerr << "[DdrReader::read] Invalid format." << std::endl;
        return 0;
    }

    if (inPhiDegrees.empty()) {
        inPhiDegrees.push_back(0.0f);
    }

    int numSpecPhi = spPhiDegrees.size();
    if (symmetryType == ddr_sdr_utility::PLANE_SYMMETRICAL) {
        numSpecPhi = spPhiDegrees.size() + (spPhiDegrees.size() - 1);
    }

    // Initialize BRDF.
    SpecularCoordinatesBrdf* brdf = new SpecularCoordinatesBrdf(inThetaDegrees.size(), inPhiDegrees.size(),
                                                                spThetaDegrees.size(), numSpecPhi,
                                                                colorModel,
                                                                numWavelengths);
    brdf->setSourceType(sourceType);

    SampleSet* ss = brdf->getSampleSet();

    copyArray(inThetaDegrees, &ss->getAngles0());
    copyArray(inPhiDegrees,   &ss->getAngles1());
    copyArray(spThetaDegrees, &ss->getAngles2());

    ss->getAngles0() = toRadians(ss->getAngles0());
    ss->getAngles1() = toRadians(ss->getAngles1());
    ss->getAngles2() = toRadians(ss->getAngles2());

    for (int i = 0; i < static_cast<int>(spPhiDegrees.size()); ++i) {
        brdf->setSpecPhi(i, toRadian(spPhiDegrees.at(i)));
    }

    // Copy symmetrical angles.
    if (symmetryType == ddr_sdr_utility::PLANE_SYMMETRICAL) {
        int numSpecPhiDegrees = static_cast<int>(spPhiDegrees.size());
        for (int i = numSpecPhiDegrees, reverseIndex = numSpecPhiDegrees - 2;
             i < brdf->getNumSpecPhi();
             ++i, --reverseIndex) {
            brdf->setSpecPhi(i, PI_F + (PI_F - brdf->getSpecPhi(reverseIndex)));
        }
    }

    std::vector<float> kbdfs;
    kbdfs.reserve(inThetaDegrees.size() * inPhiDegrees.size() * numWavelengths);

    // Read data.
    int wlIndex = 0;
    std::string dataStr;
    while (ifs >> dataStr) {
        ddr_sdr_utility::ignoreCommentLines(ifs);

        if (dataStr.empty()) {
            continue;
        }
        else if (isEqual(dataStr, "wl") ||
                 isEqual(dataStr, "bw") ||
                 isEqual(dataStr, "red") ||
                 isEqual(dataStr, "gre") ||
                 isEqual(dataStr, "blu") ||
                 isEqual(dataStr, "green") ||   // Not correct specification
                 isEqual(dataStr, "blue")) {    // Not correct specification
            if (colorModel == SPECTRAL_MODEL) {
                float wavelength;
                ifs >> wavelength;
                brdf->getSampleSet()->setWavelength(wlIndex, wavelength);
            }

            ddr_sdr_utility::ignoreCommentLines(ifs);
            // Read "kbdf" and "def" or skip "def".
            std::string kbdfStr;
            ifs >> kbdfStr;
            if (isEqual(kbdfStr, "kbdf")) {
                for (size_t i = 0; i < inThetaDegrees.size() * inPhiDegrees.size(); ++i) {
                    float kbdf;
                    ifs >> kbdf;
                    kbdfs.push_back(kbdf);
                }

                ddr_sdr_utility::ignoreCommentLines(ifs);
                // Skip "def".
                std::string defStr;
                ifs >> defStr;
            }

            int numInTheta  = static_cast<int>(inThetaDegrees.size());
            int numInPhi    = static_cast<int>(inPhiDegrees.size());
            int numSpTheta  = static_cast<int>(spThetaDegrees.size());
            int numSpPhi    = static_cast<int>(spPhiDegrees.size());

            for (int inPhIndex = 0; inPhIndex < numInPhi;   ++inPhIndex) {
            for (int inThIndex = 0; inThIndex < numInTheta; ++inThIndex) {
            for (int spPhIndex = 0; spPhIndex < numSpPhi;   ++spPhIndex) {
                ddr_sdr_utility::ignoreCommentLines(ifs);
            for (int spThIndex = 0; spThIndex < numSpTheta; ++spThIndex) {
                std::string brdfValueStr;
                ifs >> brdfValueStr;

                char* end;
                float brdfValue = static_cast<float>(std::strtod(brdfValueStr.c_str(), &end));
                if (*end != '\0') {
                    std::cerr << "[DdrReader::read] Invalid value: " << brdfValueStr << std::endl;
                    delete brdf;
                    return 0;
                }

                // Convert intensity to radiance.
                if (unitType == ddr_sdr_utility::INTENSITY_ABSOLUTE ||
                    unitType == ddr_sdr_utility::INTENSITY_RELATIVE) {
                    std::vector<float> angles;
                    angles.push_back(toRadian(inThetaDegrees.at(inThIndex)));
                    angles.push_back(toRadian(inPhiDegrees.at(inPhIndex)));
                    angles.push_back(toRadian(spThetaDegrees.at(spThIndex)));
                    angles.push_back(toRadian(spPhiDegrees.at(spPhIndex)));

                    convertCoordinateSystem<SpecularCoordinateSystem, SphericalCoordinateSystem>(
                        angles.at(0), angles.at(1), angles.at(2), angles.at(3),
                        &angles[0], &angles[1], &angles[2], &angles[3]);

                    brdfValue /= std::max(std::cos(angles.at(2)), EPSILON_F);
                    brdfValue *=  PI_F;
                }

                if (unitType == ddr_sdr_utility::LUMINANCE_ABSOLUTE ||
                    unitType == ddr_sdr_utility::INTENSITY_ABSOLUTE) {
                    if (!kbdfs.empty()) {
                        brdfValue *= kbdfs.at(inThIndex +
                                              numInTheta * inPhIndex +
                                              numInTheta * numInPhi * wlIndex);
                    }
                }

                brdfValue /= PI_F;

                Spectrum& sp = brdf->getSpectrum(inThIndex, inPhIndex, spThIndex, spPhIndex);
                sp[wlIndex] = brdfValue;

                if (symmetryType == ddr_sdr_utility::PLANE_SYMMETRICAL) {
                    int symmetryIndex = (brdf->getNumSpecPhi() - 1) - spPhIndex;
                    Spectrum& symmetrySp = brdf->getSpectrum(inThIndex, inPhIndex, spThIndex, symmetryIndex);
                    symmetrySp[wlIndex] = brdfValue;
                }

                if (ifs.fail()) {
                    std::cerr << "[DdrReader::read] Invalid format: " << brdfValue << std::endl;
                    delete brdf;
                    return 0;
                }
            }}}}

            ++wlIndex;
        }

        if (ifs.fail()) {
            std::cerr << "[DdrReader::read] Invalid format. Head of line: " << dataStr << std::endl;
            delete brdf;
            return 0;
        }
    }

    brdf->clampAngles();

    // Equalize reflectances to "kbdf"s.
    if (unitType == ddr_sdr_utility::LUMINANCE_RELATIVE ||
        unitType == ddr_sdr_utility::INTENSITY_RELATIVE) {
        if (!kbdfs.empty()) {
            Integrator integrator(PoissonDiskDistributionOnSphere::NUM_SAMPLES_ON_HEMISPHERE, true);

            int numWl = brdf->getSampleSet()->getNumWavelengths();
            int numInTh = brdf->getNumInTheta();
            int numInPh = brdf->getNumInPhi();

            for (int wlIndex = 0; wlIndex < numWl; ++wlIndex) {
                for (int inThIndex = 0; inThIndex < numInTh; ++inThIndex) {
                for (int inPhIndex = 0; inPhIndex < numInPh; ++inPhIndex) {
                    Vec3 inDir = SphericalCoordinateSystem::toXyz(brdf->getInTheta(inThIndex),
                                                                  brdf->getInPhi(inPhIndex));
                    Spectrum refSp = integrator.computeReflectance(*brdf, inDir);

                    // Edit samples with "kbdf".
                    float maxReflectance = refSp.maxCoeff();
                    for (int i2 = 0; i2 < ss->getNumAngles2(); ++i2) {
                    for (int i3 = 0; i3 < ss->getNumAngles3(); ++i3) {
                        Spectrum& sp = ss->getSpectrum(inThIndex, inPhIndex, i2, i3);

                        const float coeff = 0.999546f; // Reflectance of Lambertian using lb::Integrator.
                        sp /= maxReflectance / coeff;

                        // A reflectance equals "kbdf".
                        float kbdf = kbdfs.at(inThIndex +
                                              numInTh * inPhIndex +
                                              numInTh * numInPh * wlIndex);
                        sp *= kbdf;
                    }}
                }}
            }
        }
    }

    return brdf;
}
