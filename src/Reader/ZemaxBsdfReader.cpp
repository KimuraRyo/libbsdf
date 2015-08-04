// =================================================================== //
// Copyright (C) 2015 Kimura Ryo                                       //
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
    std::ifstream fin(fileName.c_str());
    if (fin.fail()) {
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

    ignoreCommentLines(fin);
    std::ifstream::pos_type pos = fin.tellg();

    // Read a header.
    std::string headStr;
    while (fin >> headStr) {
        ignoreCommentLines(fin);

        if (headStr.empty()) {
            continue;
        }
        else if (headStr == "Symmetry") {
            std::string typeStr;
            fin >> typeStr;

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
            fin >> typeStr;

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
            fin >> scatterStr;

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
            fin >> numInPhi;
            for (int i = 0; i < numInPhi; ++i) {
                float angle;
                fin >> angle;
                inPhiDegrees.push_back(angle);
            }
        }
        else if (headStr == "AngleOfIncidence") {
            int numInTheta;
            fin >> numInTheta;
            for (int i = 0; i < numInTheta; ++i) {
                float angle;
                fin >> angle;
                inThetaDegrees.push_back(angle);
            }
        }
        else if (headStr == "ScatterAzimuth") {
            int numSpPhi;
            fin >> numSpPhi;
            for (int i = 0; i < numSpPhi; ++i) {
                float angle;
                fin >> angle;
                spPhiDegrees.push_back(angle);
            }
        }
        else if (headStr == "ScatterRadial") {
            int numSpTheta;
            fin >> numSpTheta;
            for (int i = 0; i < numSpTheta; ++i) {
                float angle;
                fin >> angle;
                spThetaDegrees.push_back(angle);
            }
        }
        else if (headStr == "Monochrome" ||
                 headStr == "TristimulusX" ||
                 headStr == "TristimulusY" ||
                 headStr == "TristimulusZ") {
            fin.seekg(pos, std::ios_base::beg);
            break;
        }

        pos = fin.tellg();
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

    int numSpecPhi = spPhiDegrees.size();
    if (symmetryType == PLANE_SYMMETRICAL) {
        numSpecPhi = spPhiDegrees.size() + (spPhiDegrees.size() - 1);
    }

    // Initialize BRDF.
    SpecularCoordinatesBrdf* brdf = new SpecularCoordinatesBrdf(inThetaDegrees.size(), inPhiDegrees.size(),
                                                                spThetaDegrees.size(), numSpecPhi,
                                                                colorModel);

    SampleSet* ss = brdf->getSampleSet();

    copyArray(inThetaDegrees, ss->getAngles0());
    copyArray(inPhiDegrees,   ss->getAngles1());
    copyArray(spThetaDegrees, ss->getAngles2());

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
    int inThIndex = 0;
    int inPhIndex = 0;
    std::string dataStr;
    while (fin >> dataStr) {
        ignoreCommentLines(fin);

        if (dataStr.empty()) {
            continue;
        }
        else if (dataStr == "Monochrome" ||
                 dataStr == "TristimulusX") {
            wlIndex = 0;
            inThIndex = 0;
        }
        else if (dataStr == "TristimulusY") {
            wlIndex = 1;
            inThIndex = 0;
        }
        else if (dataStr == "TristimulusZ") {
            wlIndex = 2;
            inThIndex = 0;
        }
        else if (dataStr == "TIS") {
            reader_utility::ignoreLine(fin);

            for (int spPhIndex = 0; spPhIndex < static_cast<int>(spPhiDegrees.size());   ++spPhIndex) {
            for (int spThIndex = 0; spThIndex < static_cast<int>(spThetaDegrees.size()); ++spThIndex) {
                float brdfValue;
                fin >> brdfValue;

                Spectrum& sp = brdf->getSpectrum(inThIndex, inPhIndex, spThIndex, spPhIndex);
                sp[wlIndex] = brdfValue;

                if (symmetryType == PLANE_SYMMETRICAL) {
                    int symmetryIndex = (brdf->getNumSpecPhi() - 1) - spPhIndex;
                    Spectrum& symmetrySp = brdf->getSpectrum(inThIndex, inPhIndex, spThIndex, symmetryIndex);
                    symmetrySp[wlIndex] = brdfValue;
                }
            }}

            ++inThIndex;
        }

        if (fin.fail()) {
            std::cerr << "[ZemaxBsdfReader::read] Invalid format. Head of line: " << dataStr << std::endl;
            delete brdf;
            return 0;
        }
    }

    brdf->clampAngles();
    
    return brdf;
}
