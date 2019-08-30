// =================================================================== //
// Copyright (C) 2014-2019 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#include <libbsdf/Reader/SdrReader.h>

#include <fstream>

#include <libbsdf/Common/Log.h>
#include <libbsdf/Reader/DdrSdrUtility.h>

using namespace lb;

SampleSet2D* SdrReader::read(const std::string& fileName)
{
    // std::ios_base::binary is used to read line endings of CR+LF and LF.
    std::ifstream ifs(fileName.c_str(), std::ios_base::binary);
    if (ifs.fail()) {
        lbError << "[SdrReader::read] Could not open: " << fileName;
        return 0;
    }

    std::ios_base::sync_with_stdio(false);

    SourceType sourceType = UNKNOWN_SOURCE;
    ColorModel colorModel = RGB_MODEL;

    int numWavelengths = 1;
    std::vector<float> inThetaDegrees;

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
        else if (headStr == "TypeColorModel") {
            std::string typeStr;
            ifs >> typeStr;
            typeStr = reader_utility::toLower(typeStr);

            if (typeStr == "rgb") {
                colorModel = RGB_MODEL;
            }
            else if (typeStr == "spectral") {
                colorModel = SPECTRAL_MODEL;
                ifs >> numWavelengths;
            }
            else if (typeStr == "bw") {
                colorModel = MONOCHROMATIC_MODEL;
            }
            else {
                reader_utility::logNotImplementedKeyword(typeStr);
                return 0;
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
        else if (reader_utility::toLower(headStr) == "wl" ||
                 reader_utility::toLower(headStr) == "bw" ||
                 reader_utility::toLower(headStr) == "red" ||
                 reader_utility::toLower(headStr) == "gre" ||
                 reader_utility::toLower(headStr) == "blu") {
            ifs.seekg(pos, std::ios_base::beg);
            break;
        }

        pos = ifs.tellg();
    }

    if (inThetaDegrees.empty()) {
        lbError << "[SdrReader::read] Invalid format.";
        return 0;
    }

    int numInTheta = static_cast<int>(inThetaDegrees.size());

    // Initialize the array of reflectance.
    SampleSet2D* ss2 = new SampleSet2D(numInTheta, 1, colorModel, numWavelengths);

    ss2->setSourceType(sourceType);

    copyArray(inThetaDegrees, &ss2->getThetaArray());
    ss2->getThetaArray() = toRadians(ss2->getThetaArray());

    ss2->setPhi(0, 0.0);

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
                 reader_utility::toLower(dataStr) == "blu") {
            if (colorModel == SPECTRAL_MODEL) {
                float wavelength;
                ifs >> wavelength;
                ss2->setWavelength(wlIndex, wavelength);
            }

            ddr_sdr_utility::ignoreCommentLines(ifs);
            // Skip "def"
            std::string defStr;
            ifs >> defStr;

            for (int inThIndex = 0; inThIndex < numInTheta; ++inThIndex) {
                std::string srValueStr;
                ifs >> srValueStr;

                char* end;
                float srValue = static_cast<float>(std::strtod(srValueStr.c_str(), &end));
                if (*end != '\0') {
                    lbError << "[SdrReader::read] Invalid value: " << srValueStr;
                    delete ss2;
                    return 0;
                }

                ss2->getSpectrum(inThIndex)[wlIndex] = srValue;

                if (ifs.fail()) {
                    lbError << "[SdrReader::read] Invalid format: " << srValue;
                    delete ss2;
                    return 0;
                }
            }

            ++wlIndex;
        }

        if (ifs.fail()) {
            lbError << "[SdrReader::read] Invalid format. Head of line: " << dataStr;
            delete ss2;
            return 0;
        }
    }

    ss2->clampAngles();

    return ss2;
}
