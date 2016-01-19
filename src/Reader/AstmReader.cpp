// =================================================================== //
// Copyright (C) 2014-2015 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#include <libbsdf/Reader/AstmReader.h>

#include <fstream>
#include <iostream>
#include <set>
#include <sstream>

#include <libbsdf/Brdf/Processor.h>
#include <libbsdf/Brdf/RandomSampleSet.h>

using namespace lb;

SphericalCoordinatesBrdf* AstmReader::read(const std::string& fileName)
{
    std::ifstream fin(fileName.c_str(), std::ifstream::binary);
    if (fin.fail()) {
        std::cerr << "[AstmReader::read] Could not open: " << fileName << std::endl;
        return 0;
    }

    std::ios_base::sync_with_stdio(false);

    ColorModel colorModel = SPECTRAL_MODEL;
    std::vector<float> wavelengths;

    // Read a header.
    std::string headStr;
    while (fin >> headStr) {
        if (headStr.empty()) {
            continue;
        }
        else if (headStr == "NUM_POINTS") {
            int numPoints;
            fin >> numPoints;
            std::cout << "[AstmReader::read] NUM_POINTS: " << numPoints << std::endl;
        }
        else if (headStr == "VARS") {
            fin.ignore(1);

            std::string varStr;
            std::getline(fin, varStr);
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
                return 0;
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

    std::set<float> inThetaAngles;
    std::set<float> inPhiAngles;
    std::set<float> outThetaAngles;
    std::set<float> outPhiAngles;

    RandomSampleSet rss;
    RandomSampleSet::SampleMap& samples = rss.getSampleMap();

    // Read data.
    std::string dataStr;
    while (std::getline(fin, dataStr)) {
        if (dataStr.empty() || dataStr.at(0) == '\r') continue;

        RandomSampleSet::AngleList angles;
        Spectrum values(wavelengths.size());

        std::stringstream stream(dataStr);
        std::string token;
        int count = 0;
        while (std::getline(stream, token, ',')) {
            if (token.at(token.size() - 1) == '\r') {
                token.erase(token.size() - 1);
            }

            std::stringstream stream(token);
            float val;
            stream >> val;

            if (count <= 3) {
                if (val < 0.0f) {
                    val += 2.0f * PI_F;
                }

                val = std::min(val, SphericalCoordinateSystem::MAX_ANGLE3);
                angles.push_back(val);
            }
            else {
                values[count - 4] = std::max(val, 0.0f);
            }

            ++count;
        }

        inThetaAngles.insert(angles.at(0));
        inPhiAngles.insert(angles.at(1));
        outThetaAngles.insert(angles.at(2));
        outPhiAngles.insert(angles.at(3));

        auto it = samples.find(angles);
        if (it != samples.end()) {
            std::cout
                << "[AstmReader::read] Already defined: "
                << angles.at(0) << ", " << angles.at(1) << ", " << angles.at(2) << ", " << angles.at(3)
                << std::endl;
            continue;
        }
        
        samples[angles] = values;
    }

    // Modify data for the isotropic BRDF with the incoming azimuthal angle of non-zero radian.
    if (inPhiAngles.size() == 1 && *inPhiAngles.begin() != 0.0f) {
        // Rotate outgoing azimuthal angles using the incoming azimuthal angle.
        RandomSampleSet::AngleList angles;
        for (auto it = outPhiAngles.begin(); it != outPhiAngles.end(); ++it) {
            float outPhi = *it - *inPhiAngles.begin();
            if (outPhi < 0.0f) {
                outPhi += 2.0f * PI_F;
            }

            angles.push_back(outPhi);
        }
        outPhiAngles.clear();
        std::copy(angles.begin(), angles.end(), std::inserter(outPhiAngles, outPhiAngles.begin()));

        // Rotate sample points using the incoming azimuthal angle.
        RandomSampleSet::SampleMap modifiedSamples;
        for (auto it = samples.begin(); it != samples.end(); ++it) {
            RandomSampleSet::AngleList angles = it->first;
            float outPhi = angles.at(3) - *inPhiAngles.begin();
            if (outPhi < 0.0f) {
                outPhi += 2.0f * PI_F;
            }
            angles.at(1) = 0.0f;
            angles.at(3) = outPhi;

            modifiedSamples[angles] = it->second;
        }
        samples.clear();
        std::copy(modifiedSamples.begin(), modifiedSamples.end(), std::inserter(samples, samples.begin()));

        inPhiAngles.clear();
        inPhiAngles.insert(0.0f);
    }

    outPhiAngles.insert(0.0f);
    outPhiAngles.insert(SphericalCoordinateSystem::MAX_ANGLE3);

    int numInTheta  = inThetaAngles.size();
    int numInPhi    = inPhiAngles.size();
    int numOutTheta = outThetaAngles.size();
    int numOutPhi   = outPhiAngles.size();
    SphericalCoordinatesBrdf* brdf = new SphericalCoordinatesBrdf(numInTheta, numInPhi,
                                                                  numOutTheta, numOutPhi,
                                                                  colorModel, wavelengths.size());
    SampleSet* ss = brdf->getSampleSet();

    // Set angles.
    copyArray(inThetaAngles,  &ss->getAngles0());
    copyArray(inPhiAngles,    &ss->getAngles1());
    copyArray(outThetaAngles, &ss->getAngles2());
    copyArray(outPhiAngles,   &ss->getAngles3());

    // Set wavelengths.
    for (int i = 0; i < static_cast<int>(wavelengths.size()); ++i) {
        ss->setWavelength(i, wavelengths.at(i));
    }

    rss.setupBrdf(brdf);

    std::cout << "[AstmReader::read] One side of the plane of incidence: " << ss->isOneSide() << std::endl;

    if (ss->isOneSide()) {
        SphericalCoordinatesBrdf* filledBrdf = fillSymmetricBrdf(brdf);
        delete brdf;
        brdf = filledBrdf;
    }

    std::cout << "[AstmReader::read] The number of sample points: " << samples.size() << std::endl;

    brdf->clampAngles();

    return brdf;
}
