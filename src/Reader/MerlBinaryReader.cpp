// =================================================================== //
// Copyright (C) 2014-2019 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#include <libbsdf/Reader/MerlBinaryReader.h>

#include <fstream>

#include <libbsdf/Common/Log.h>
#include <libbsdf/Common/Vector.h>

using namespace lb;

HalfDifferenceCoordinatesBrdf* MerlBinaryReader::read(const std::string& fileName)
{
    std::ifstream ifs(fileName.c_str(), std::ios_base::binary);
    if (ifs.fail()) {
        lbError << "[MerlBinaryReader::read] Could not open: " << fileName;
        return 0;
    }

    std::ios_base::sync_with_stdio(false);

    const int numHalfTheta = 90;
    const int numDiffTheta = 90;
    const int numDiffPhi = 360;
    const int numSamples = numHalfTheta * numDiffTheta * numDiffPhi / 2;

    // Read a header.
    int dims[3];
    ifs.read(reinterpret_cast<char*>(dims), sizeof(int) * 3);

    int numSamplesInFile = dims[0] * dims[1] * dims[2];
    if (numSamplesInFile != numSamples) {
        lbError << "[MerlBinaryReader::read] Dimensions do not match: " << numSamplesInFile << ", " << numSamples;
        return 0;
    }

    // Read data.
    double* samples = new double[numSamples * 3];
    ifs.read(reinterpret_cast<char*>(samples), sizeof(double) * numSamples * 3);
    if (ifs.fail()) {
        lbError << "[MerlBinaryReader::read] Invalid format.";
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

    const Vec3f rgbScaleCoeff(1.0f / 1500.0f, 1.15f / 1500.0f, 1.66f / 1500.0f);

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

        Vec3f rgb(static_cast<float>(samples[sampleIndex]),
                  static_cast<float>(samples[sampleIndex + numSamples]),
                  static_cast<float>(samples[sampleIndex + numSamples * 2]));
        rgb = rgb.cwiseMax(0.0f);

        rgb = rgb.cwiseProduct(rgbScaleCoeff);
        brdf->setSpectrum(halfThIndex, 0, diffThIndex, diffPhIndex, rgb);
    }}}

    delete[] samples;

    brdf->clampAngles();
    brdf->setSourceType(MEASURED_SOURCE);

    return brdf;
}
