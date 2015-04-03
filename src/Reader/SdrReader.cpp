// =================================================================== //
// Copyright (C) 2014-2015 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#include <libbsdf/Reader/SdrReader.h>

#include <fstream>
#include <iostream>
#include <limits>

#include <libbsdf/Reader/DdrSdrUtility.h>

using namespace lb;

SampleSet2D* SdrReader::read(const std::string& fileName)
{
    std::ifstream fin(fileName.c_str());
    if (fin.fail()) {
        std::cerr << "[SdrReader::read] Could not open: " << fileName << std::endl;
        return 0;
    }

    SampleSet2D* ss2 = 0;

    ColorModel::Type colorModel = ColorModel::RGB;

    int numWavelengths = 1;
    std::vector<float> inThetaDegrees;
    int wlIndex = 0;

    ddr_sdr_utility::ignoreCommentLines(fin);

    std::string headStr;
    while (fin >> headStr) {
        ddr_sdr_utility::ignoreCommentLines(fin);

        if (headStr.empty()) {
            continue;
        }
        else if (headStr == "Source") {
            reader_utility::ignoreLine(fin);
        }
        else if (headStr == "TypeColorModel") {
            std::string typeStr;
            fin >> typeStr;
            typeStr = reader_utility::toLower(typeStr);

            if (typeStr == "rgb") {
                colorModel = ColorModel::RGB;
            }
            else if (typeStr == "spectral") {
                colorModel = ColorModel::SPECTRAL;
                fin >> numWavelengths;
            }
            else if (typeStr == "bw") {
                colorModel = ColorModel::MONOCHROME;
            }
        }
        else if (headStr == "sigma") {
            int numInTheta;
            fin >> numInTheta;
            for (int i = 0; i < numInTheta; ++i) {
                float angle;
                fin >> angle;
                inThetaDegrees.push_back(angle);
            }
        }
        else if (reader_utility::toLower(headStr) == "wl") {
            if (!ss2) {
                if (inThetaDegrees.empty()) {
                    inThetaDegrees.push_back(0.0);
                }

                ss2 = new SampleSet2D(inThetaDegrees.size(), 1, colorModel, numWavelengths);

                copy(inThetaDegrees, ss2->getThetaArray());
                ss2->getThetaArray() = toRadians(ss2->getThetaArray());

                ss2->setPhi(0, 0.0);
            }

            float wavelength;
            fin >> wavelength;
            ss2->setWavelength(wlIndex, wavelength);

            ddr_sdr_utility::ignoreCommentLines(fin);
            // Skip "def"
            std::string defStr;
            fin >> defStr;

            for (int inThIndex = 0; inThIndex < static_cast<int>(inThetaDegrees.size()); ++inThIndex) {
                float srValue;
                fin >> srValue;

                ss2->getSpectrum(inThIndex)[wlIndex] = srValue;

                if (fin.fail()) {
                    std::cerr << "[SdrReader::read] Invalid format: " << srValue << std::endl;
                    delete ss2;
                    return 0;
                }
            }

            ++wlIndex;
        }
        else {
            //logNotImplementedKeyword(headStr);
        }

        if (fin.fail()) {
            std::cerr << "[SdrReader::read] Invalid format. Head of line: " << headStr << std::endl;
            delete ss2;
            return 0;
        }
    }

    ss2->clampAngles();
    
    return ss2;
}
