// =================================================================== //
// Copyright (C) 2014-2023 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#include <libbsdf/Reader/AstmReader.h>

#include <fstream>
#include <iterator>
#include <set>
#include <sstream>

#include <libbsdf/Brdf/Processor.h>
#include <libbsdf/Brdf/SphericalCoordinatesRandomSampleSet.h>

using namespace lb;

SphericalCoordinatesBrdf* AstmReader::read(const std::string& fileName)
{
    // std::ios_base::binary is used to read line endings of CR+LF and LF.
    std::ifstream ifs(fileName.c_str(), std::ios_base::binary);
    if (ifs.fail()) {
        lbError << "[AstmReader::read] Could not open: " << fileName;
        return nullptr;
    }

    std::ios_base::sync_with_stdio(false);

    ColorModel colorModel = SPECTRAL_MODEL;
    std::vector<float> wavelengths;

    // Read a header.
    std::string headStr;
    while (ifs >> headStr) {
        if (headStr.empty()) {
            continue;
        }
        else if (headStr == "NUM_POINTS") {
            int numPoints;
            ifs >> numPoints;
            lbInfo << "[AstmReader::read] NUM_POINTS: " << numPoints;
        }
        else if (headStr == "VARS") {
            ifs.ignore(1);

            std::string varStr;
            std::getline(ifs, varStr);
            std::stringstream stream(varStr);
            std::string token;
            std::vector<std::string> varNames;
            while (getline(stream, token, ',')) {
                if (token.at(token.size() - 1) == '\r') {
                    token.erase(token.size() - 1);
                }

                varNames.push_back(token);
            }

            // Read the names of color components.
            if (varNames.size() < 5) {
                // No color components.
                return nullptr;
            }
            else if (varNames.size() == 7 &&
                     varNames.at(4) == "R" &&
                     varNames.at(5) == "G" &&
                     varNames.at(6) == "B") {
                wavelengths.push_back(0.0f);
                wavelengths.push_back(0.0f);
                wavelengths.push_back(0.0f);
                colorModel = RGB_MODEL;
            }
            else if (varNames.size() == 7 &&
                     varNames.at(4) == "X" &&
                     varNames.at(5) == "Y" &&
                     varNames.at(6) == "Z") {
                wavelengths.push_back(0.0f);
                wavelengths.push_back(0.0f);
                wavelengths.push_back(0.0f);
                colorModel = XYZ_MODEL;
            }
            else {
                for (int i = 4; i < static_cast<int>(varNames.size()); ++i) {
                    std::string name = varNames.at(i);
                    bool spectral = (name.size() >= 3 &&
                                     name.substr(name.size() - 2, 2) == "nm");
                    if (spectral) {
                        std::stringstream nameStream(name.substr(0, name.size() - 2));
                        float wl;
                        nameStream >> wl;

                        if (nameStream.fail()) {
                            wavelengths.push_back(0.0f);
                        }
                        else {
                            wavelengths.push_back(wl);
                        }
                    }
                    else {
                        wavelengths.push_back(0.0f);
                    }
                }
            }

            break;
        }
    }

    if (colorModel == SPECTRAL_MODEL) {
        if (wavelengths.size() == 1 && wavelengths.at(0) == 0.0f) {
            colorModel = MONOCHROMATIC_MODEL;
        }
        else if (wavelengths.size() == 3 && wavelengths.at(0) == 0.0f) {
            colorModel = RGB_MODEL;
        }
    }

    std::set<double> inThetaAngles;
    std::set<double> inPhiAngles;
    std::set<double> outThetaAngles;
    std::set<double> outPhiAngles;

    SphericalCoordinatesRandomSampleSet rss;
    SphericalCoordinatesRandomSampleSet::SampleMap& samples = rss.getSampleMap();

    // Read data.
    std::string dataStr;
    while (std::getline(ifs, dataStr)) {
        if (dataStr.empty() || dataStr.at(0) == '\r') continue;

        SphericalCoordinatesRandomSampleSet::AngleList angles;
        Spectrum values(wavelengths.size());

        std::stringstream stream(dataStr);
        std::string token;
        int index = 0;
        while (std::getline(stream, token, ',')) {
            if (token.at(token.size() - 1) == '\r') {
                token.erase(token.size() - 1);
            }

            double val = std::atof(token.c_str());

            if (index <= 3) {
                if (val < 0) {
                    val += TAU_D;
                }

                val = std::min(val, SphericalCoordinateSystem::MAX_ANGLE3);
                angles.push_back(val);
            }
            else {
                int valIndex = index - 4;
                if (wavelengths.size() == valIndex) {
                    lbError << "[AstmReader::read] Data format is invalid.";
                    return nullptr;
                }

                values[valIndex] = static_cast<float>(std::max(val, 0.0));
            }

            ++index;
        }

        if (wavelengths.size() < index - 4) {
            lbError << "[AstmReader::read] Data format is invalid.";
            return nullptr;
        }

        inThetaAngles.insert(angles.at(0));
        inPhiAngles.insert(angles.at(1));
        outThetaAngles.insert(angles.at(2));
        outPhiAngles.insert(angles.at(3));

        auto it = samples.find(angles);
        if (it != samples.end()) {
            lbWarn
                << "[AstmReader::read] Already defined: "
                << angles.at(0) << ", " << angles.at(1) << ", " << angles.at(2) << ", " << angles.at(3);
            continue;
        }

        samples[angles] = values;
    }

    // Modify data for the isotropic BRDF with the incoming azimuthal angle of non-zero radian.
    if (inPhiAngles.size() == 1 && *inPhiAngles.begin() != 0) {
        // Rotate outgoing azimuthal angles using the incoming azimuthal angle.
        SphericalCoordinatesRandomSampleSet::AngleList rotatedAngles;
        for (auto it = outPhiAngles.begin(); it != outPhiAngles.end(); ++it) {
            double outPhi = *it - *inPhiAngles.begin();
            if (outPhi < 0) {
                outPhi += TAU_D;
            }

            rotatedAngles.push_back(outPhi);
        }
        outPhiAngles.clear();
        std::copy(rotatedAngles.begin(), rotatedAngles.end(), std::inserter(outPhiAngles, outPhiAngles.begin()));

        // Rotate sample points using the incoming azimuthal angle.
        SphericalCoordinatesRandomSampleSet::SampleMap modifiedSamples;
        for (auto it = samples.begin(); it != samples.end(); ++it) {
            SphericalCoordinatesRandomSampleSet::AngleList angles = it->first;
            double outPhi = angles.at(3) - *inPhiAngles.begin();
            if (outPhi < 0) {
                outPhi += TAU_D;
            }
            angles.at(1) = 0;
            angles.at(3) = outPhi;

            modifiedSamples[angles] = it->second;
        }
        samples.clear();
        std::copy(modifiedSamples.begin(), modifiedSamples.end(), std::inserter(samples, samples.begin()));

        inPhiAngles.clear();
        inPhiAngles.insert(0);
    }

    outPhiAngles.insert(0);
    outPhiAngles.insert(SphericalCoordinateSystem::MAX_ANGLE3);

    SphericalCoordinatesBrdf* brdf = new SphericalCoordinatesBrdf(static_cast<int>(inThetaAngles.size()),
                                                                  static_cast<int>(inPhiAngles.size()),
                                                                  static_cast<int>(outThetaAngles.size()),
                                                                  static_cast<int>(outPhiAngles.size()),
                                                                  colorModel,
                                                                  static_cast<int>(wavelengths.size()));
    SampleSet* ss = brdf->getSampleSet();

    // Set angles.
    array_util::copy(inThetaAngles,  &ss->getAngles0());
    array_util::copy(inPhiAngles,    &ss->getAngles1());
    array_util::copy(outThetaAngles, &ss->getAngles2());
    array_util::copy(outPhiAngles,   &ss->getAngles3());

    // Set wavelengths.
    for (int i = 0; i < static_cast<int>(wavelengths.size()); ++i) {
        ss->setWavelength(i, wavelengths.at(i));
    }

    rss.setupBrdf(brdf);

    lbInfo << "[AstmReader::read] One side of the plane of incidence: " << ss->isOneSide();

    if (ss->isOneSide()) {
        SphericalCoordinatesBrdf* filledBrdf = fillSymmetricBrdf(brdf);
        delete brdf;
        brdf = filledBrdf;
    }

    lbInfo << "[AstmReader::read] The number of sample points: " << samples.size();

    brdf->clampAngles();
    brdf->setSourceType(MEASURED_SOURCE);

    return brdf;
}
