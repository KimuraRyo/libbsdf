// =================================================================== //
// Copyright (C) 2014-2017 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#include <libbsdf/Reader/DdrReader.h>

#include <fstream>

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

    // Read a header.
    std::string headStr;
    while (ifs >> headStr) {
        ddr_sdr_utility::ignoreCommentLines(ifs);

        if (headStr.empty()) {
            continue;
        }
        else if (headStr == "Source") {
            std::string typeStr;
            ifs >> typeStr;

            if (typeStr == "Measured") {
                sourceType = MEASURED_SOURCE;
            }
            else if (typeStr == "Generated") {
                sourceType = GENERATED_SOURCE;
            }
            else if (typeStr == "Edited") {
                sourceType = EDITED_SOURCE;
            }
            else if (typeStr == "Morphed") {
                sourceType = UNKNOWN_SOURCE;
            }
            else {
                reader_utility::logNotImplementedKeyword(typeStr);
            }
        }
        else if (headStr == "TypeSym") {
            std::string typeStr;
            ifs >> typeStr;
            typeStr = reader_utility::toLower(typeStr);

            if (typeStr == reader_utility::toLower("AxiSymmetrical")) {
                symmetryType = ddr_sdr_utility::AXI_SYMMETRICAL;
            }
            else if (typeStr == reader_utility::toLower("DirSymmetrical")) {
                symmetryType = ddr_sdr_utility::DIR_SYMMETRICAL;
            }
            else if (typeStr == reader_utility::toLower("PlaneSymmetrical")) {
                symmetryType = ddr_sdr_utility::PLANE_SYMMETRICAL;
            }
            else if (typeStr == reader_utility::toLower("ASymmetrical")) {
                std::ifstream::pos_type asymmetricalPos = ifs.tellg();
                ifs >> typeStr;
                typeStr = reader_utility::toLower(typeStr);

                if (typeStr == reader_utility::toLower("4D")) {                
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
        else if (headStr == "TypeColorModel") {
            std::string typeStr;
            ifs >> typeStr;
            typeStr = reader_utility::toLower(typeStr);

            if (typeStr == "rgb") {
                colorModel = RGB_MODEL;
                numWavelengths = 3;
            }
            else if (typeStr == "spectral") {
                colorModel = SPECTRAL_MODEL;
                ifs >> numWavelengths;
            }
            else if (typeStr == "bw") {
                colorModel = MONOCHROMATIC_MODEL;
                numWavelengths = 1;
            }
            else {
                reader_utility::logNotImplementedKeyword(typeStr);
                return 0;
            }
        }
        else if (headStr == "TypeData") {
            std::string typeStr;
            ifs >> typeStr;

            if (typeStr == "Luminance") {
                ifs >> typeStr;

                if (typeStr == "Absolute") {
                    unitType = ddr_sdr_utility::LUMINANCE_ABSOLUTE;
                }
                else if (typeStr == "Relative") {
                    unitType = ddr_sdr_utility::LUMINANCE_RELATIVE;
                }
                else {
                    reader_utility::ignoreLine(ifs);
                }
            }
            else {
                reader_utility::ignoreLine(ifs);
            }
        }
        else if (headStr == "psi") {
            int numInPhi;
            ifs >> numInPhi;
            for (int i = 0; i < numInPhi; ++i) {
                float angle;
                ifs >> angle;
                inPhiDegrees.push_back(angle);
            }
        }
        else if (headStr == "sigma") {
            int numInTheta;
            ifs >> numInTheta;
            for (int i = 0; i < numInTheta; ++i) {
                float angle;
                ifs >> angle;
                inThetaDegrees.push_back(angle);
            }
        }
        else if (headStr == "sigmaT") {
            reader_utility::logNotImplementedKeyword(headStr);
            reader_utility::ignoreLine(ifs);
        }
        else if (headStr == "phi") {
            int numSpecPhi;
            ifs >> numSpecPhi;
            for (int i = 0; i < numSpecPhi; ++i) {
                float angle;
                ifs >> angle;
                spPhiDegrees.push_back(angle);
            }
        }
        else if (headStr == "theta") {
            int numSpecTheta;
            ifs >> numSpecTheta;
            for (int i = 0; i < numSpecTheta; ++i) {
                float angle;
                ifs >> angle;
                spThetaDegrees.push_back(angle);
            }
        }
        else if (reader_utility::toLower(headStr) == "wl" ||
                 reader_utility::toLower(headStr) == "bw" ||
                 reader_utility::toLower(headStr) == "red" ||
                 reader_utility::toLower(headStr) == "gre" ||
                 reader_utility::toLower(headStr) == "blu" ||
                 reader_utility::toLower(headStr) == "green" || // Not correct specification
                 reader_utility::toLower(headStr) == "blue") {  // Not correct specification
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

    // Read data.
    int wlIndex = 0;
    std::string dataStr;
    while (ifs >> dataStr) {
        ddr_sdr_utility::ignoreCommentLines(ifs);

        if (dataStr.empty()) {
            continue;
        }
        else if (reader_utility::toLower(dataStr) == "wl" ||
                 reader_utility::toLower(dataStr) == "bw" ||
                 reader_utility::toLower(dataStr) == "red" ||
                 reader_utility::toLower(dataStr) == "gre" ||
                 reader_utility::toLower(dataStr) == "blu" ||
                 reader_utility::toLower(dataStr) == "green" || // Not correct specification
                 reader_utility::toLower(dataStr) == "blue") {  // Not correct specification
            if (colorModel == SPECTRAL_MODEL) {
                float wavelength;
                ifs >> wavelength;
                brdf->getSampleSet()->setWavelength(wlIndex, wavelength);
            }

            ddr_sdr_utility::ignoreCommentLines(ifs);
            // Read "kbdf" and "def" or skip "def".
            std::string kbdfStr;
            ifs >> kbdfStr;
            std::vector<float> kbdfs;
            if (reader_utility::toLower(kbdfStr) == "kbdf") {
                for (int i = 0;
                     i < static_cast<int>(inThetaDegrees.size() * inPhiDegrees.size());
                     ++i) {
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
                float brdfValue = static_cast<float>(std::atof(brdfValueStr.c_str()));

                if (unitType == ddr_sdr_utility::LUMINANCE_ABSOLUTE ||
                    unitType == ddr_sdr_utility::LUMINANCE_RELATIVE) {
                    brdfValue /= PI_F;
                }

                if (!kbdfs.empty()) {
                    brdfValue *= kbdfs.at(inThIndex + numInTheta * inPhIndex);
                }

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

    return brdf;
}
