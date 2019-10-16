// =================================================================== //
// Copyright (C) 2019 Kimura Ryo                                       //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#include <iostream>
#include <memory>

#include <libbsdf/Brdf/Processor.h>

#include <libbsdf/Reader/AstmReader.h>
#include <libbsdf/Reader/DdrReader.h>
#include <libbsdf/Reader/ReaderUtility.h>

#include <libbsdf/Writer/DdrWriter.h>

#include <ArgumentParser.h>
#include <Utility.h>

using namespace lb;

/*
 * BRDF inserter
 */

const std::string APP_NAME("lbinsert");
const std::string APP_VERSION("1.0.1");

void showHelp()
{
    using std::cout;
    using std::endl;

    cout << "Usage: lbinsert [options ...] partial_file angle base_file" << endl;
    cout << endl;
    cout << "lbinsert inserts a BRDF with an incoming azimuthal angle into a base BRDF." << endl;
    cout << "Two files must have the same type of color." << endl;
    cout << endl;
    cout << "Positional Arguments:" << endl;
    cout << "  partial_file     Name of a BRDF file inserted into a base BRDF." << endl;
    cout << "                   This BRDF must have an incoming azimuthal angle." << endl;
    cout << "                   \".ddr\" or \".astm\" is acceptable as a suffix." << endl;
    cout << "  angle            Incoming azimuthal angle in degrees (range: [0, 360])" << endl;
    cout << "  base_file        Name of a base BRDF file." << endl;
    cout << "                   \".ddr\" is acceptable as a suffix." << endl;
    cout << endl;
    cout << "Options:" << endl;
    cout << "  -h, --help       show this help message and exit" << endl;
    cout << "  -v, --version    show program's version number and exit" << endl;
}

int main(int argc, char** argv)
{
    Log::setNotificationLevel(Log::Level::WARN_MSG);

    ArgumentParser ap(argc, argv);

    if (ap.read("-h") || ap.read("--help") || ap.getTokens().empty()) {
        showHelp();
        return 0;
    }

    if (ap.read("-v") || ap.read("--version")) {
        app_utility::showAppVersion(APP_NAME, APP_VERSION);
        return 0;
    }

    if (!ap.validateNumTokens(3)) return 1;

    std::string partFileName    = ap.getTokens().at(0);
    std::string angleStr        = ap.getTokens().at(1);
    std::string baseFileName    = ap.getTokens().at(2);

    FileType partFileType = reader_utility::classifyFile(partFileName);
    FileType baseFileType = reader_utility::classifyFile(baseFileName);

    // Load a base BRDF file.
    std::unique_ptr<SpecularCoordinatesBrdf> baseBrdf;
    if (baseFileType == INTEGRA_DDR_FILE) {
        baseBrdf.reset(DdrReader::read(baseFileName));

        if (!baseBrdf) {
            std::cerr << "Failed to load: " << baseFileName << std::endl;
            return 1;
        }
    }
    else {
        std::cerr << "Unsupported file type: " << baseFileType << std::endl;
        return 1;
    }

    // Load a partial BRDF file.
    std::unique_ptr<SpecularCoordinatesBrdf> partBrdf;
    switch (partFileType) {
        case ASTM_FILE: {
            std::unique_ptr<SphericalCoordinatesBrdf> brdf(AstmReader::read(partFileName));
            if (!brdf) {
                std::cerr << "Failed to load: " << partFileName << std::endl;
                return 1;
            }

            int numPolarAngles      = baseBrdf->getNumSpecTheta();
            int numAzimuthalAngles  = baseBrdf->getNumSpecPhi();
            partBrdf.reset(new SpecularCoordinatesBrdf(*brdf, numPolarAngles, numAzimuthalAngles));
            break;
        }
        case INTEGRA_DDR_FILE:
            partBrdf.reset(DdrReader::read(partFileName));
            break;
        default:
            std::cerr << "Unsupported file type: " << partFileType << std::endl;
            return 1;
    }

    if (!partBrdf) {
        std::cerr << "Failed to load: " << partFileName << std::endl;
        return 1;
    }

    // Validate color models.
    const SampleSet* baseSs = baseBrdf->getSampleSet();
    SampleSet* partSs       = partBrdf->getSampleSet();
    if (!hasSameColor(*baseSs, *partSs)) {
        if (partSs->getColorModel()      == SPECTRAL_MODEL &&
            partSs->getNumWavelengths()  == 1) {
            // If two monochromatic data sets have different color models, partial data is adjusted.
            partSs->setColorModel(MONOCHROMATIC_MODEL);
            partSs->setWavelength(0, 0.0f);
        }
        else {
            std::cerr << "Color models or wavelengths do not match." << std::endl;
            return 1;
        }
    }

    // Read an incoming azimuthal angle.
    char* end;
    double inAzimuthalDegree = std::strtod(angleStr.c_str(), &end);
    if (*end != '\0') {
        std::cerr << "Invalid value: " << angleStr << std::endl;
        return 1;
    }
    inAzimuthalDegree = app_utility::clampParameter("angle", inAzimuthalDegree, 0.0, 360.0);

    // Insert the partial BRDF to the base BRDF.
    float inAzimuthalAngle = static_cast<float>(lb::toRadian(inAzimuthalDegree));
    std::unique_ptr<SpecularCoordinatesBrdf> brdf(insertBrdfAlongInPhi(*baseBrdf, *partBrdf, inAzimuthalAngle));

    if (!brdf) {
        std::cerr << "Failed to insert a BRDF." << std::endl;
    }

    // Save a DDR file.
    std::string comments = app_utility::createComments(argc, argv, APP_NAME, APP_VERSION);
    if (DdrWriter::write(baseFileName, *brdf, comments)) {
        std::cout << "Saved: " << baseFileName << std::endl;
    }

    return 0;
}
