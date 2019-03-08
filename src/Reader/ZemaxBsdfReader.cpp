// =================================================================== //
// Copyright (C) 2015-2019 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#include <libbsdf/Reader/ZemaxBsdfReader.h>

#include <fstream>
#include <iostream>
#include <set>

using namespace lb;

SpecularCoordinatesBrdf* ZemaxBsdfReader::read(const std::string& fileName, DataType* dataType)
{
    // std::ios_base::binary is used to read line endings of CR+LF and LF.
    std::ifstream ifs(fileName.c_str(), std::ios_base::binary);
    if (ifs.fail()) {
        std::cerr << "[ZemaxBsdfReader::read] Could not open: " << fileName << std::endl;
        return 0;
    }

    std::ios_base::sync_with_stdio(false);

    SymmetryType symmetryType = UNKNOWN_SYMMETRY;
    ColorModel colorModel = UNKNOWN_MODEL;

    std::vector<float> inThetaDegrees;
    std::vector<float> inPhiDegrees;
    std::vector<float> spThetaDegrees;
    std::vector<float> spPhiDegrees;

    int numChannels = 1;

    ignoreCommentLines(ifs);
    std::ifstream::pos_type pos = ifs.tellg();

    // Read a header.
    std::string headStr;
    while (ifs >> headStr) {
        ignoreCommentLines(ifs);

        if (headStr.empty()) {
            continue;
        }
        else if (headStr == "Symmetry") {
            std::string typeStr;
            ifs >> typeStr;

            if (typeStr == "PlaneSymmetrical") {
                symmetryType = PLANE_SYMMETRICAL;
            }
            else if (typeStr == "Asymmetrical") {
                symmetryType = ASYMMETRICAL;
            }
            else if (typeStr == "ASymmetrical4D") {
                symmetryType = ASYMMETRICAL_4D;
            }
            else {
                reader_utility::logNotImplementedKeyword(typeStr);
                return 0;
            }
        }
        else if (headStr == "SpectralContent") {
            std::string typeStr;
            ifs >> typeStr;

            if (typeStr == "Monochrome") {
                colorModel = MONOCHROMATIC_MODEL;
                numChannels = 1;
            }
            else if (typeStr == "XYZ") {
                colorModel = XYZ_MODEL;
                numChannels = 3;
            }
            else {
                reader_utility::logNotImplementedKeyword(typeStr);
                return 0;
            }
        }
        else if (headStr == "ScatterType") {
            std::string scatterStr;
            ifs >> scatterStr;

            if (scatterStr == "BRDF") {
                *dataType = BRDF_DATA;
            }
            else if (scatterStr == "BTDF") {
                *dataType = BTDF_DATA;
            }
            else {
                reader_utility::logNotImplementedKeyword(scatterStr);
                return 0;
            }
        }
        else if (headStr == "SampleRotation") {
            int numInPhi;
            ifs >> numInPhi;
            for (int i = 0; i < numInPhi; ++i) {
                float angle;
                ifs >> angle;
                inPhiDegrees.push_back(angle);
            }
        }
        else if (headStr == "AngleOfIncidence") {
            int numInTheta;
            ifs >> numInTheta;
            for (int i = 0; i < numInTheta; ++i) {
                float angle;
                ifs >> angle;
                inThetaDegrees.push_back(angle);
            }
        }
        else if (headStr == "ScatterAzimuth") {
            int numSpPhi;
            ifs >> numSpPhi;
            for (int i = 0; i < numSpPhi; ++i) {
                float angle;
                ifs >> angle;
                spPhiDegrees.push_back(angle);
            }
        }
        else if (headStr == "ScatterRadial") {
            int numSpTheta;
            ifs >> numSpTheta;
            for (int i = 0; i < numSpTheta; ++i) {
                float angle;
                ifs >> angle;
                spThetaDegrees.push_back(angle);
            }
        }
        else if (headStr == "Monochrome" ||
                 headStr == "TristimulusX" ||
                 headStr == "TristimulusY" ||
                 headStr == "TristimulusZ") {
            ifs.seekg(pos, std::ios_base::beg);
            break;
        }

        pos = ifs.tellg();
    }

    if (inThetaDegrees.empty() ||
        spThetaDegrees.empty() ||
        spPhiDegrees.empty()) {
        std::cerr << "[ZemaxBsdfReader::read] Invalid format." << std::endl;
        return 0;
    }

    if (inPhiDegrees.empty()) {
        inPhiDegrees.push_back(0.0f);
    }

    size_t numSpecPhi = spPhiDegrees.size();
    if (symmetryType == PLANE_SYMMETRICAL) {
        numSpecPhi = spPhiDegrees.size() + (spPhiDegrees.size() - 1);
    }

    // Initialize BRDF.
    SpecularCoordinatesBrdf* brdf = new SpecularCoordinatesBrdf(static_cast<int>(inThetaDegrees.size()),
                                                                static_cast<int>(inPhiDegrees.size()),
                                                                static_cast<int>(spThetaDegrees.size()),
                                                                static_cast<int>(numSpecPhi),
                                                                colorModel);
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
    if (symmetryType == PLANE_SYMMETRICAL) {
        int numSpecPhiDegrees = static_cast<int>(spPhiDegrees.size());
        for (int i = numSpecPhiDegrees, reverseIndex = numSpecPhiDegrees - 2;
             i < brdf->getNumSpecPhi();
             ++i, --reverseIndex) {
            brdf->setSpecPhi(i, PI_F + (PI_F - brdf->getSpecPhi(reverseIndex)));
        }
    }

    // Read data.
    int wlIndex = 0;
    int cntTis = 0;
    std::string dataStr;
    while (ifs >> dataStr) {
        ignoreCommentLines(ifs);

        if (dataStr.empty()) {
            continue;
        }
        else if (dataStr == "Monochrome" ||
                 dataStr == "TristimulusX") {
            wlIndex = 0;
            cntTis = 0;
        }
        else if (dataStr == "TristimulusY") {
            wlIndex = 1;
            cntTis = 0;
        }
        else if (dataStr == "TristimulusZ") {
            wlIndex = 2;
            cntTis = 0;
        }
        else if (dataStr == "TIS") {
            reader_utility::ignoreLine(ifs);

            int inPhIndex = cntTis / static_cast<int>(inThetaDegrees.size());
            int inThIndex = cntTis - static_cast<int>(inThetaDegrees.size()) * inPhIndex;
            
            for (int spPhIndex = 0; spPhIndex < static_cast<int>(spPhiDegrees.size());   ++spPhIndex) {
            for (int spThIndex = 0; spThIndex < static_cast<int>(spThetaDegrees.size()); ++spThIndex) {
                std::string brdfValueStr;
                ifs >> brdfValueStr;
                float brdfValue = static_cast<float>(std::atof(brdfValueStr.c_str()));

                Spectrum& sp = brdf->getSpectrum(inThIndex, inPhIndex, spThIndex, spPhIndex);
                sp[wlIndex] = brdfValue;

                if (symmetryType == PLANE_SYMMETRICAL) {
                    int symmetryIndex = (brdf->getNumSpecPhi() - 1) - spPhIndex;
                    Spectrum& symmetrySp = brdf->getSpectrum(inThIndex, inPhIndex, spThIndex, symmetryIndex);
                    symmetrySp[wlIndex] = brdfValue;
                }
            }}

            ++cntTis;
        }

        if (ifs.fail()) {
            std::cerr << "[ZemaxBsdfReader::read] Invalid format. Head of line: " << dataStr << std::endl;
            delete brdf;
            return 0;
        }
    }

    brdf->clampAngles();
    brdf->setSourceType(MEASURED_SOURCE);
    
    return brdf;
}
