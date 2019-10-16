// =================================================================== //
// Copyright (C) 2019 Kimura Ryo                                       //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#include <iostream>
#include <memory>

#include <libbsdf/Reader/AstmReader.h>
#include <libbsdf/Reader/DdrReader.h>
#include <libbsdf/Reader/ZemaxBsdfReader.h>

#include <libbsdf/Writer/DdrWriter.h>

#include <ArgumentParser.h>
#include <Utility.h>

using namespace lb;

/*
 * BRDF/BTDF converter
 */

const std::string APP_NAME("lbconv");
const std::string APP_VERSION("1.0.1");

void showHelp()
{
    using std::cout;
    using std::endl;

    cout << "Usage: lbconv [options ...] in_file out_file" << endl;
    cout << endl;
    cout << "lbconv converts a BRDF/BTDF file to an Integra BRDF/BTDF file." << endl;
    cout << endl;
    cout << "Positional Arguments:" << endl;
    cout << "  in_file      Name of an input BRDF/BTDF file." << endl;
    cout << "               Valid formats:" << endl;
    cout << "                   Integra Diffuse Distribution (\".ddr, .ddt\")" << endl;
    cout << "                   Zemax BSDF (\".bsdf\")" << endl;
    cout << "                   ASTM E1392-96(2002) (\".astm\")" << endl;
    cout << "  out_file     Name of an output BRDF/BTDF file." << endl;
    cout << "               \".ddr\" for BRDF or \".ddt\" for BTDF is acceptable as a suffix." << endl;
    cout << "               If an appropriate suffix is not obtained, \".ddr\" or \".ddt\" is appended." << endl;
    cout << endl;
    cout << "Options:" << endl;
    cout << "  -h, --help       show this help message and exit" << endl;
    cout << "  -v, --version    show program's version number and exit" << endl;
    cout << "  -scatterType     set either BRDF or BTDF for the input ASTM file (default: BRDF)" << endl;
    cout << "  -arrangement     arrange BRDF/BTDF with extrapolation and conservation of energy" << endl;
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

    DataType dataType = BRDF_DATA;
    std::string dataTypeStr;
    if (ap.read("-scatterType", &dataTypeStr)) {
        if (dataTypeStr == "BRDF") {
            dataType = BRDF_DATA;
        }
        else if (dataTypeStr == "BTDF") {
            dataType = BTDF_DATA;
        }
        else {
            std::cerr << "Invalid scatter type: " << dataTypeStr << std::endl;
            return 1;
        }
    }

    bool arranged = false;
    if (ap.read("-arrangement")) {
        arranged = true;
    }

    if (!ap.validateNumTokens(2)) return 1;

    std::string inFileName  = ap.getTokens().at(0);
    std::string outFileName = ap.getTokens().at(1);

    // Load a BRDF/BTDF file.
    FileType inFileType = reader_utility::classifyFile(inFileName);
    std::unique_ptr<Brdf> inBrdf;
    switch (inFileType) {
        case ASTM_FILE:
            inBrdf.reset(AstmReader::read(inFileName));
            break;
        case INTEGRA_DDR_FILE:
            inBrdf.reset(DdrReader::read(inFileName));
            dataType = BRDF_DATA;
            break;
        case INTEGRA_DDT_FILE:
            inBrdf.reset(DdrReader::read(inFileName));
            dataType = BTDF_DATA;
            break;
        case lb::ZEMAX_FILE:
            inBrdf.reset(ZemaxBsdfReader::read(inFileName, &dataType));
            break;
        default:
            std::cerr << "Unsupported file type: " << inFileType << std::endl;
            return 1;
    }

    if (!inBrdf) {
        std::cerr << "Failed to load: " << inFileName << std::endl;
        return 1;
    }

    // Convert the BRDF/BTDF.
    std::unique_ptr<SpecularCoordinatesBrdf> outBrdf(DdrWriter::convert(*inBrdf));
    if (arranged) {
        outBrdf.reset(DdrWriter::arrange(*outBrdf, dataType));
    }

    // Fix the output filename.
    if (dataType == BRDF_DATA && !reader_utility::hasSuffix(outFileName, ".ddr")) {
        outFileName += ".ddr";
    }
    else if (dataType == BTDF_DATA && !reader_utility::hasSuffix(outFileName, ".ddt")) {
        outFileName += ".ddt";
    }

    // Save a DDR/DDT file.
    std::string comments = app_utility::createComments(argc, argv, APP_NAME, APP_VERSION);
    if (DdrWriter::write(outFileName, *outBrdf, comments)) {
        std::cout << "Saved: " << outFileName << std::endl;
    }

    return 0;
}
