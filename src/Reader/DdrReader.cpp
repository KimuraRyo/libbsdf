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

    ddr_sdr_utility::SymmetryType symmetryType  = ddr_sdr_utility::PLANE_SYMMETRICAL;
    ddr_sdr_utility::UnitType     unitType      = ddr_sdr_utility::LUMINANCE_ABSOLUTE;

    ColorModel colorModel = RGB_MODEL;

    int numWavelengths = 1;

    std::vector<float> inThetaDegrees;
    std::vector<float> inPhiDegrees;
    std::vector<float> spThetaDegrees;
    std::vector<float> spPhiDegrees;

    ddr_sdr_utility::ignoreCommentLines(fin);
    std::ifstream::pos_type pos = fin.tellg();

    // Read a header.
    std::string headStr;
    while (fin >> headStr) {
        ddr_sdr_utility::ignoreCommentLines(fin);

        if (headStr.empty()) {
            continue;
        }
        else if (headStr == "Source") {
            reader_utility::ignoreLine(fin);
        }
        else if (headStr == "TypeSym") {
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
                return 0;
            }
        }
        else if (headStr == "TypeColorModel") {
            std::string typeStr;
            fin >> typeStr;
            typeStr = reader_utility::toLower(typeStr);

            if (typeStr == "rgb") {
                colorModel = RGB_MODEL;
                numWavelengths = 3;
            }
            else if (typeStr == "spectral") {
                colorModel = SPECTRAL_MODEL;
                fin >> numWavelengths;
            }
            else if (typeStr == "bw") {
                colorModel = MONOCHROMATIC_MODEL;
                numWavelengths = 1;
            }
            else {
                reader_utility::logNotImplementedKeyword(typeStr);
                return 0;
            }
        }
        else if (headStr == "TypeData") {
            std::string typeStr;
            fin >> typeStr;

            if (typeStr == "Luminance Absolute") {
                unitType = ddr_sdr_utility::LUMINANCE_ABSOLUTE;
            }
            else if (typeStr == "Luminance Relative") {
                unitType = ddr_sdr_utility::LUMINANCE_RELATIVE;
            }
            else {
                reader_utility::ignoreLine(fin);
            }
        }
        else if (headStr == "psi") {
            int numInPhi;
            fin >> numInPhi;
            for (int i = 0; i < numInPhi; ++i) {
                float angle;
                fin >> angle;
                inPhiDegrees.push_back(angle);
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
        else if (headStr == "sigmaT") {
            reader_utility::logNotImplementedKeyword(headStr);
            reader_utility::ignoreLine(fin);
        }
        else if (headStr == "phi") {
            int numSpecPhi;
            fin >> numSpecPhi;
            for (int i = 0; i < numSpecPhi; ++i) {
                float angle;
                fin >> angle;
                spPhiDegrees.push_back(angle);
            }
        }
        else if (headStr == "theta") {
            int numSpecTheta;
            fin >> numSpecTheta;
            for (int i = 0; i < numSpecTheta; ++i) {
                float angle;
                fin >> angle;
                spThetaDegrees.push_back(angle);
            }
        }
        else if (reader_utility::toLower(headStr) == "wl" ||
                 reader_utility::toLower(headStr) == "bw" ||
                 reader_utility::toLower(headStr) == "red" ||
                 reader_utility::toLower(headStr) == "gre" ||
                 reader_utility::toLower(headStr) == "blu") {
            fin.seekg(pos, std::ios_base::beg);
            break;
        }

        pos = fin.tellg();
    }

    if (inThetaDegrees.empty() ||
        spThetaDegrees.empty() ||
        spPhiDegrees.empty()) {
        std::cerr << "[DdrReader::read] Invalid format." << std::endl;
        return 0;
    }

    if (inPhiDegrees.empty()) {
        inPhiDegrees.push_back(0.0f);
    }

    int numSpecPhi = spPhiDegrees.size();
    if (symmetryType == ddr_sdr_utility::PLANE_SYMMETRICAL) {
        numSpecPhi = spPhiDegrees.size() + (spPhiDegrees.size() - 1);
    }

    // Initialize BRDF.
    SpecularCoordinatesBrdf* brdf = new SpecularCoordinatesBrdf(inThetaDegrees.size(), inPhiDegrees.size(),
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

    // Copy symmetrical angles.
    if (symmetryType == ddr_sdr_utility::PLANE_SYMMETRICAL) {
        int numSpecPhiDegrees = static_cast<int>(spPhiDegrees.size());
        for (int i = numSpecPhiDegrees, reverseIndex = numSpecPhiDegrees - 2;
             i < brdf->getNumSpecPhi();
             ++i, --reverseIndex) {
            brdf->setSpecPhi(i, PI_F + (PI_F - brdf->getSpecPhi(reverseIndex)));
        }
    }

    // Read data.
    int wlIndex = 0;
    std::string dataStr;
    while (fin >> dataStr) {
        ddr_sdr_utility::ignoreCommentLines(fin);

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
                fin >> wavelength;
                brdf->getSampleSet()->setWavelength(wlIndex, wavelength);
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

                if (unitType == ddr_sdr_utility::LUMINANCE_ABSOLUTE ||
                    unitType == ddr_sdr_utility::LUMINANCE_RELATIVE) {
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

        if (fin.fail()) {
            std::cerr << "[DdrReader::read] Invalid format. Head of line: " << dataStr << std::endl;
            delete brdf;
            return 0;
        }
    }

    brdf->clampAngles();

    return brdf;
}
