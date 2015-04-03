// =================================================================== //
// Copyright (C) 2014-2015 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#include <libbsdf/Reader/DdrReader.h>

#include <libbsdf/Reader/DdrSdrUtility.h>

using namespace lb;

SpecularCoordinatesBrdf* DdrReader::read(const std::string& fileName)
{
    std::ifstream fin(fileName.c_str());
    if (fin.fail()) {
        std::cerr << "[DdrReader::read] Could not open: " << fileName << std::endl;
        return 0;
    }

    std::ios_base::sync_with_stdio(false);

    SpecularCoordinatesBrdf* brdf = 0;

    ddr_sdr_utility::SymmetryType symmetryType  = ddr_sdr_utility::PLANE_SYMMETRICAL;
    ddr_sdr_utility::DataType     dataType      = ddr_sdr_utility::LUMINANCE_ABSOLUTE;

    ColorModel::Type colorModel = ColorModel::RGB;

    int numWavelengths = 1;

    std::vector<float> inThetaDegrees;
    std::vector<float> inPhiDegrees;
    std::vector<float> spThetaDegrees;
    std::vector<float> spPhiDegrees;

    int wlIndex = 0;

    ddr_sdr_utility::ignoreCommentLines(fin);

    std::string token;
    while (fin >> token) {
        ddr_sdr_utility::ignoreCommentLines(fin);

        if (token.empty()) {
            continue;
        }
        else if (token == "Source") {
            reader_utility::ignoreLine(fin);
        }
        else if (token == "TypeSym") {
            std::string typeStr;
            fin >> typeStr;

            if (typeStr == "AxiSymmetrical") {
                symmetryType = ddr_sdr_utility::AXI_SYMMETRICAL;
            }
            else if (typeStr == "DirSymmetrical") {
                symmetryType = ddr_sdr_utility::DIR_SYMMETRICAL;
            }
            else if (typeStr == "PlaneSymmetrical") {
                symmetryType = ddr_sdr_utility::PLANE_SYMMETRICAL;
            }
            else if (typeStr == "ASymmetrical") {
                symmetryType = ddr_sdr_utility::ASYMMETRICAL;
            }
            else if (typeStr == "ASymmetrical 4D") {
                symmetryType = ddr_sdr_utility::ASYMMETRICAL_4D;
            }
            else {
                reader_utility::logNotImplementedKeyword(typeStr);
            }
        }
        else if (token == "TypeColorModel") {
            std::string typeStr;
            fin >> typeStr;
            typeStr = reader_utility::toLower(typeStr);

            if (typeStr == "rgb") {
                colorModel = ColorModel::RGB;
                numWavelengths = 3;
            }
            else if (typeStr == "spectral") {
                colorModel = ColorModel::SPECTRAL;
                fin >> numWavelengths;
            }
            else if (typeStr == "bw") {
                colorModel = ColorModel::MONOCHROME;
                numWavelengths = 1;
            }
        }
        else if (token == "TypeData") {
            std::string typeStr;
            fin >> typeStr;

            if (typeStr == "Luminance Absolute") {
                dataType = ddr_sdr_utility::LUMINANCE_ABSOLUTE;
            }
            else if (typeStr == "Luminance Relative") {
                dataType = ddr_sdr_utility::LUMINANCE_RELATIVE;
            }
            else {
                reader_utility::ignoreLine(fin);
            }
        }
        else if (token == "psi") {
            int numInPhi;
            fin >> numInPhi;
            for (int i = 0; i < numInPhi; ++i) {
                float angle;
                fin >> angle;
                inPhiDegrees.push_back(angle);
            }
        }
        else if (token == "sigma") {
            int numInTheta;
            fin >> numInTheta;
            for (int i = 0; i < numInTheta; ++i) {
                float angle;
                fin >> angle;
                inThetaDegrees.push_back(angle);
            }
        }
        else if (token == "sigmaT") {
            reader_utility::logNotImplementedKeyword(token);
            reader_utility::ignoreLine(fin);
        }
        else if (token == "phi") {
            int numSpecPhi;
            fin >> numSpecPhi;
            for (int i = 0; i < numSpecPhi; ++i) {
                float angle;
                fin >> angle;
                spPhiDegrees.push_back(angle);
            }
        }
        else if (token == "theta") {
            int numSpecTheta;
            fin >> numSpecTheta;
            for (int i = 0; i < numSpecTheta; ++i) {
                float angle;
                fin >> angle;
                spThetaDegrees.push_back(angle);
            }
        }
        else if (reader_utility::toLower(token) == "wl" ||
                 reader_utility::toLower(token) == "bw" ||
                 reader_utility::toLower(token) == "red" ||
                 reader_utility::toLower(token) == "gre" ||
                 reader_utility::toLower(token) == "blu") {
            if (!brdf) {
                if (inThetaDegrees.empty()) {
                    inThetaDegrees.push_back(0.0);
                }

                if (inPhiDegrees.empty()) {
                    inPhiDegrees.push_back(0.0);
                }

                if (spThetaDegrees.empty()) {
                    spThetaDegrees.push_back(0.0);
                }

                if (spPhiDegrees.empty()) {
                    spPhiDegrees.push_back(0.0);
                }

                int numSpecPhi = spPhiDegrees.size();
                if (symmetryType == ddr_sdr_utility::PLANE_SYMMETRICAL) {
                    numSpecPhi = spPhiDegrees.size() + (spPhiDegrees.size() - 1);
                }

                if (spThetaDegrees.size() <= 1 || numSpecPhi <= 1) return 0;

                brdf = new SpecularCoordinatesBrdf(inThetaDegrees.size(), inPhiDegrees.size(),
                                                   spThetaDegrees.size(), numSpecPhi,
                                                   colorModel,
                                                   numWavelengths);

                SampleSet* ss = brdf->getSampleSet();

                copy(inThetaDegrees, ss->getAngles0());
                copy(inPhiDegrees,   ss->getAngles1());
                copy(spThetaDegrees, ss->getAngles2());

                ss->getAngles0() = toRadians(ss->getAngles0());
                ss->getAngles1() = toRadians(ss->getAngles1());
                ss->getAngles2() = toRadians(ss->getAngles2());

                for (int i = 0; i < static_cast<int>(spPhiDegrees.size()); ++i) {
                    brdf->setSpecPhi(i, toRadian(spPhiDegrees.at(i)));
                }
                
                // Copy symmetry samples.
                if (symmetryType == ddr_sdr_utility::PLANE_SYMMETRICAL) {
                    int numSpecPhiDegrees = static_cast<int>(spPhiDegrees.size());
                    for (int i = numSpecPhiDegrees, reverseIndex = numSpecPhiDegrees - 2;
                         i < brdf->getNumSpecPhi();
                         ++i, --reverseIndex) {
                        brdf->setSpecPhi(i, PI_F + (PI_F - brdf->getSpecPhi(reverseIndex)));
                    }
                }
            }

            if (colorModel == ColorModel::SPECTRAL) {
                float wavelength;
                fin >> wavelength;
                brdf->getSampleSet()->setWavelength(wlIndex, wavelength);
            }
            else if (colorModel == ColorModel::RGB) {
                brdf->getSampleSet()->setWavelength(wlIndex, 0.0f);
            }
            else {
                brdf->getSampleSet()->setWavelength(0, 0.0f);
            }

            ddr_sdr_utility::ignoreCommentLines(fin);
            std::string kbdfStr;
            fin >> kbdfStr;
            std::vector<float> kbdfs;
            if (reader_utility::toLower(kbdfStr) == "kbdf") {
                for (int i = 0; i < static_cast<int>(inThetaDegrees.size()); ++i) {
                    float kbdf;
                    fin >> kbdf;
                    kbdfs.push_back(kbdf);
                }
            }

            ddr_sdr_utility::ignoreCommentLines(fin);
            // Skip "def"
            std::string defStr;
            fin >> defStr;

            for (int inPhIndex = 0; inPhIndex < static_cast<int>(inPhiDegrees.size());   ++inPhIndex) {
            for (int inThIndex = 0; inThIndex < static_cast<int>(inThetaDegrees.size()); ++inThIndex) {
            for (int spPhIndex = 0; spPhIndex < static_cast<int>(spPhiDegrees.size());   ++spPhIndex) {
                ddr_sdr_utility::ignoreCommentLines(fin);
            for (int spThIndex = 0; spThIndex < static_cast<int>(spThetaDegrees.size()); ++spThIndex) {
                float brdfValue;
                fin >> brdfValue;

                if (dataType == ddr_sdr_utility::LUMINANCE_ABSOLUTE ||
                    dataType == ddr_sdr_utility::LUMINANCE_RELATIVE) {
                    brdfValue /= PI_F;
                }

                if (!kbdfs.empty()) {
                    brdfValue *= kbdfs.at(inThIndex);
                }

                Spectrum& sp = brdf->getSpectrum(inThIndex, inPhIndex, spThIndex, spPhIndex);
                sp[wlIndex] = brdfValue;

                if (symmetryType == ddr_sdr_utility::PLANE_SYMMETRICAL) {
                    int symmetryIndex = (brdf->getNumSpecPhi() - 1) - spPhIndex;
                    Spectrum& symmetrySp = brdf->getSpectrum(inThIndex, inPhIndex, spThIndex, symmetryIndex);
                    symmetrySp[wlIndex] = brdfValue;
                }

                if (fin.fail()) {
                    std::cerr << "[DdrReader::read] Invalid format: " << brdfValue << std::endl;
                    delete brdf;
                    return 0;
                }
            }}}}

            ++wlIndex;
        }
        else {
            //logNotImplementedKeyword(token);
        }

        if (fin.fail()) {
            std::cerr << "[DdrReader::read] Invalid format. Head of line: " << token << std::endl;
            delete brdf;
            return 0;
        }
    }

    brdf->clampAngles();

    return brdf;
}
