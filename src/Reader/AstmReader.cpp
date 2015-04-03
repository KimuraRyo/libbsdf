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

using namespace lb;

SphericalCoordinatesBrdf* AstmReader::read(const std::string& fileName)
{
    std::ifstream fin(fileName.c_str(), std::ifstream::binary);
    if (fin.fail()) {
        std::cerr << "[AstmReader::read] Could not open: " << fileName << std::endl;
        return 0;
    }

    std::ios_base::sync_with_stdio(false);

    ColorModel::Type colorModel = ColorModel::SPECTRAL;
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
                colorModel = ColorModel::RGB;
            }
            else if (varNames.size() == 7 &&
                     varNames.at(4) == "X" &&
                     varNames.at(5) == "Y" &&
                     varNames.at(6) == "Z") {
                wavelengths.push_back(0.0f);
                wavelengths.push_back(0.0f);
                wavelengths.push_back(0.0f);
                colorModel = ColorModel::XYZ;
            }
            else {
                for (int i = 4; i < static_cast<int>(varNames.size()); ++i) {
                    std::string name = varNames.at(i);
                    bool isWavelength = (name.size() >= 3 && name.substr(name.size() - 2, 2) == "nm");
                    if (isWavelength) {
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

    if (colorModel == ColorModel::SPECTRAL) {
        if (wavelengths.size() == 1 && wavelengths.at(0) == 0.0f) {
            colorModel = ColorModel::MONOCHROME;
        }
        else if (wavelengths.size() == 3 && wavelengths.at(0) == 0.0f) {
            colorModel = ColorModel::RGB;
        }
    }

    std::set<float> inThetaAngles;
    std::set<float> inPhiAngles;
    std::set<float> outThetaAngles;
    std::set<float> outPhiAngles;

    SampleMap samples;

    // Read data.
    std::string valStr;
    while (std::getline(fin, valStr)) {
        if (valStr.empty() || valStr.at(0) == '\r') continue;

        AngleList angles;
        Spectrum values(wavelengths.size());

        std::stringstream stream(valStr);
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

            count++;
        }

        inThetaAngles.insert(angles.at(0));
        inPhiAngles.insert(angles.at(1));
        outThetaAngles.insert(angles.at(2));
        outPhiAngles.insert(angles.at(3));

        auto it = samples.find(angles);
        if (it != samples.end()) {
            std::cout
                << "[AstmReader::read] Already inserted: "
                << angles.at(0) << ", " << angles.at(1) << ", " << angles.at(2) << ", " << angles.at(3) << std::endl;
            continue;
        }
        
        samples[angles] = values;
    }

    // Modify data for an isotropic BRDF.
    if (inPhiAngles.size() == 1 && *inPhiAngles.begin() != 0.0f) {
        // Modify outgoing azimuthal angles.
        {
            AngleList angles;
            for (auto it = outPhiAngles.begin(); it != outPhiAngles.end(); ++it) {
                float outPhi = *it - *inPhiAngles.begin();
                if (outPhi < 0.0f) {
                    outPhi += 2.0f * PI_F;
                }

                angles.push_back(outPhi);
            }

            outPhiAngles.clear();
            std::copy(angles.begin(), angles.end(), std::inserter(outPhiAngles, outPhiAngles.begin()));
        }

        // Modify sample points.
        {
            SampleMap modifiedSamples;
            for (auto it = samples.begin(); it != samples.end(); ++it) {
                AngleList angles = it->first;
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
        }

        inPhiAngles.clear();
        inPhiAngles.insert(0.0f);
    }

    outPhiAngles.insert(0.0f);
    outPhiAngles.insert(SphericalCoordinateSystem::MAX_ANGLE3);

    int numInTheta  = inThetaAngles.size();
    int numInPhi    = inPhiAngles.size();
    int numOutTheta = outThetaAngles.size();
    int numOutPhi   = outPhiAngles.size();
    SphericalCoordinatesBrdf* brdf = new SphericalCoordinatesBrdf(numInTheta, numInPhi, numOutTheta, numOutPhi,
                                                                  colorModel, wavelengths.size());

    SampleSet* ss = brdf->getSampleSet();

    // Set angles.
    copy(inThetaAngles,  ss->getAngles0());
    copy(inPhiAngles,    ss->getAngles1());
    copy(outThetaAngles, ss->getAngles2());
    copy(outPhiAngles,   ss->getAngles3());

    // Set wavelengths.
    for (int i = 0; i < static_cast<int>(wavelengths.size()); ++i) {
        ss->setWavelength(i, wavelengths.at(i));
    }

    bool containOutPhi_0_PI = false;
    bool containOutPhi_PI_2PI = false;

    for (int inThIndex  = 0; inThIndex  < brdf->getNumInTheta();  ++inThIndex)  {
    for (int inPhIndex  = 0; inPhIndex  < brdf->getNumInPhi();    ++inPhIndex)  {
    for (int outThIndex = 0; outThIndex < brdf->getNumOutTheta(); ++outThIndex) {
        AngleList angles;
        float outPhi;
        SampleMap::iterator it;
        #pragma omp parallel for private(angles, outPhi, it)
        for (int outPhIndex = 0; outPhIndex < brdf->getNumOutPhi(); ++outPhIndex) {
            angles.resize(4);
            angles.at(0) = brdf->getInTheta(inThIndex);
            angles.at(1) = brdf->getInPhi(inPhIndex);
            angles.at(2) = brdf->getOutTheta(outThIndex);
            angles.at(3) = brdf->getOutPhi(outPhIndex);

            // Check the omission of a plane symmetry BRDF.
            outPhi = angles.at(3);
            if (!containOutPhi_0_PI && outPhi > 0.0f && outPhi < PI_F) {
                containOutPhi_0_PI = true;
            }
            if (!containOutPhi_PI_2PI && outPhi > PI_F && outPhi < 2.0f * PI_F) {
                containOutPhi_PI_2PI = true;
            }

            it = samples.find(angles);
            if (it != samples.end()) {
                brdf->setSpectrum(inThIndex, inPhIndex, outThIndex, outPhIndex,
                                  it->second);
            }
            else {
                brdf->setSpectrum(inThIndex, inPhIndex, outThIndex, outPhIndex,
                                  findNearestSample(samples, angles));
            }
        }
    }}}

    std::cout
        << "[AstmReader::read] The omission of a plane symmetry BRDF: "
        << (!containOutPhi_0_PI || !containOutPhi_PI_2PI) << std::endl;

    if (!containOutPhi_0_PI || !containOutPhi_PI_2PI) {
        SphericalCoordinatesBrdf* filledBrdf = fillSymmetricBrdf(brdf);
        delete brdf;
        brdf = filledBrdf;
    }

    std::cout << "[AstmReader::read] The number of sample points: " << samples.size() << std::endl;

    brdf->clampAngles();

    return brdf;
}

SphericalCoordinatesBrdf* AstmReader::fillSymmetricBrdf(SphericalCoordinatesBrdf* brdf)
{
    AngleList filledAngles;

    for (int i = 0; i < brdf->getNumOutPhi(); ++i) {
        float outPhi = brdf->getOutPhi(i);
        bool isOmittedAngle = (outPhi != 0.0f &&
                               !isEqual(outPhi, PI_F) &&
                               !isEqual(outPhi, 2.0f * PI_F));
        if (isOmittedAngle) {
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
    {
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
    }

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

        int origIndex;
        // Find the corresponding index.
        for (origIndex = 0; origIndex < brdf->getNumOutPhi(); ++origIndex) {
            float origOutPhi = brdf->getOutPhi(origIndex);
            bool isEqualOutPhi = (origOutPhi == outPhi ||
                                  isEqual(origOutPhi, SphericalCoordinateSystem::MAX_ANGLE3 - outPhi));
            if (isEqualOutPhi) break;
        }

        Spectrum& sp = brdf->getSpectrum(inThIndex, inPhIndex, outThIndex, origIndex);
        filledBrdf->setSpectrum(inThIndex, inPhIndex, outThIndex, outPhIndex, sp);
    }}}}

    return filledBrdf;
}

const Spectrum& AstmReader::findNearestSample(const SampleMap& samples, const AngleList& sampleAngles)
{
    const Spectrum* sp = 0;

    float minError = std::numeric_limits<float>::max();
    for (auto it = samples.begin(); it != samples.end(); ++it) {
        const AngleList& angles = it->first;

        if (angles.at(0) != sampleAngles.at(0) ||
            angles.at(1) != sampleAngles.at(1)) continue;

        Vec3 sampleOutDir = SphericalCoordinateSystem::toXyz(sampleAngles.at(2), sampleAngles.at(3));
        Vec3 outDir = SphericalCoordinateSystem::toXyz(angles.at(2), angles.at(3));

        float error = (sampleOutDir - outDir).norm();
        if (error == 0.0f) {
            return it->second;
        }
        else if (error < minError) {
            minError = error;
            sp = &it->second;
        }
    }

    return *sp;
}
