// =================================================================== //
// Copyright (C) 2014-2015 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#include <libbsdf/Reader/MerlBinaryReader.h>

#include <fstream>
#include <iostream>

#include <libbsdf/Common/Vector.h>

using namespace lb;

HalfDifferenceCoordinatesBrdf* MerlBinaryReader::read(const std::string& fileName)
{
    std::ifstream fin(fileName.c_str(), std::ifstream::binary);
    if (fin.fail()) {
        std::cerr << "[MerlBinaryReader::read] Could not open: " << fileName << std::endl;
        return 0;
    }

    const int numHalfTheta = 90;
    const int numDiffTheta = 90;
    const int numDiffPhi = 360;
    const int numSamples = numHalfTheta * numDiffTheta * numDiffPhi / 2;

    int dims[3];
    fin.read(reinterpret_cast<char*>(dims), sizeof(int) * 3);

    int numSamplesInFile = dims[0] * dims[1] * dims[2];
    if (numSamplesInFile != numSamples) {
        std::cerr
            << "[MerlBinaryReader::read] Dimensions don't match: " << numSamplesInFile << ", " << numSamples
            << std::endl;
        return 0;
    }

    double* samples = new double[numSamples * 3];
    fin.read(reinterpret_cast<char*>(samples), sizeof(double) * numSamples * 3);
    if (fin.fail()) {
        std::cerr << "[MerlBinaryReader::read] Invalid format" << std::endl;
        delete[] samples;
    }

    HalfDifferenceCoordinatesBrdf* brdf = new HalfDifferenceCoordinatesBrdf(numHalfTheta + 1, 1,
                                                                            numDiffTheta + 1, numDiffPhi + 1,
                                                                            RGB_MODEL, 3, true);

    // Set the angles of a non-linear mapping.
    for (int i = 0; i < brdf->getNumHalfTheta(); ++i) {
        float halfThetaDegree = static_cast<float>(i * i) / numHalfTheta;
        brdf->setHalfTheta(i, toRadian(halfThetaDegree));
    }

    const Vec3 rgbScaleCoeff(1.0f / 1500.0f, 1.15f / 1500.0f, 1.66f / 1500.0f);

    for (int halfThIndex = 0; halfThIndex < brdf->getNumHalfTheta(); ++halfThIndex) {
    for (int diffThIndex = 0; diffThIndex < brdf->getNumDiffTheta(); ++diffThIndex) {
    for (int diffPhIndex = 0; diffPhIndex < brdf->getNumDiffPhi();   ++diffPhIndex) {
        // Compute the sample index including undefined samples.
        int sampleHalfThIndex = std::min(halfThIndex, numHalfTheta - 1);
        int sampleDiffThIndex = std::min(diffThIndex, numDiffTheta - 1);
        int sampleDiffPhIndex = diffPhIndex % (numDiffPhi / 2);

        int sampleIndex = sampleDiffPhIndex
                        + numDiffPhi / 2 * sampleDiffThIndex
                        + numDiffPhi / 2 * numDiffTheta * sampleHalfThIndex;

        Vec3 rgb(static_cast<float>(samples[sampleIndex]),
                 static_cast<float>(samples[sampleIndex + numSamples]),
                 static_cast<float>(samples[sampleIndex + numSamples * 2]));
        rgb = rgb.cwiseMax(0.0f);

        rgb = rgb.cwiseProduct(rgbScaleCoeff);
        brdf->setSpectrum(halfThIndex, 0, diffThIndex, diffPhIndex, rgb.asVector3f());
    }}}

    delete[] samples;

    brdf->clampAngles();

    return brdf;
}
