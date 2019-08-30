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
#include <libbsdf/Brdf/SphericalCoordinatesBrdf.h>
#include <libbsdf/Brdf/SpecularCoordinatesBrdf.h>

#include <libbsdf/Common/Log.h>
#include <libbsdf/Common/Version.h>

#include <libbsdf/Reader/AstmReader.h>
#include <libbsdf/Reader/DdrReader.h>
#include <libbsdf/Reader/ReaderUtility.h>
#include <libbsdf/Reader/ZemaxBsdfReader.h>

#include <libbsdf/Writer/DdrWriter.h>

#include <ArgumentParser.h>

using namespace lb;

/*
 * BRDF/BTDF converter
 */

int main(int argc, char** argv)
{
    Log::setNotificationLevel(Log::Level::ERROR_MSG);

    ArgumentParser ap(argc, argv);

    using std::cout;
    using std::cerr;
    using std::endl;

    if (ap.read("-h") || ap.read("--help") || ap.getTokens().empty()) {
        cout << "Usage: lbconv [options ...] in_file out_file" << endl;
        cout << endl;
        cout << "lbconv converts a BRDF/BTDF file to an Integra BRDF/BTDF file." << endl;
        cout << endl;
        cout << "Positional Arguments:" << endl;
        cout << "  in_file  Name of an input BRDF/BTDF file." << endl;
        cout << "           Valid formats:" << endl;
        cout << "               Integra Diffuse Distribution (\".ddr, .ddt\")" << endl;
        cout << "               Zemax BSDF (\".bsdf\")" << endl;
        cout << "               ASTM E1392-96(2002) (\".astm\")" << endl;
        cout << "  out_file Name of an output BRDF/BTDF file." << endl;
        cout << "           \".ddr\" for BRDF or \".ddt\" for BTDF is acceptable as a suffix." << endl;
        cout << "           If an appropriate suffix is not obtained, \".ddr\" or \".ddt\" is appended." << endl;
        cout << endl;
        cout << "Options:" << endl;
        cout << "  -h, --help       show this help message and exit" << endl;
        cout << "  -v, --version    show program's version number and exit" << endl;
        cout << "  -scatterType     set either BRDF or BTDF for the input ASTM file (default: BRDF)" << endl;
        cout << "  -arrangement     arrange BRDF/BTDF with extrapolate and conservation of energy." << endl;
        return 0;
    }

    const std::string version("1.0.0");

    if (ap.read("-v") || ap.read("--version")) {
        cout << "Version: lbconv " << version << " (libbsdf-" << getVersion() << ")" << endl;
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
            cerr << "Invalid scatter type: " << dataTypeStr << endl;
            return 1;
        }
    }

    bool arranged = false;
    if (ap.read("-arrangement")) {
        arranged = true;
    }

    if (!ap.validate(2)) return 1;

    std::string inFileName  = ap.getTokens().at(0);
    std::string outFileName = ap.getTokens().at(1);

    FileType inFileType = reader_utility::classifyFile(inFileName);

    std::unique_ptr<Brdf> inBrdf;

    // Load a BRDF/BTDF file.
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
            cerr << "Invalid file type: " << inFileName << endl;
            return 1;
    }

    if (!inBrdf) {
        cerr << "Failed to load: " << inFileName << endl;
        return 1;
    }

    if (dataType == BRDF_DATA && !reader_utility::hasSuffix(outFileName, ".ddr")) {
        outFileName += ".ddr";
    }
    else if (dataType == BTDF_DATA && !reader_utility::hasSuffix(outFileName, ".ddt")) {
        outFileName += ".ddt";
    }

    std::string comments("Software: lbconv-" + version);
    comments += "\n;; Arguments:";
    for (int i = 1; i < argc; ++i) {
        comments += " " + std::string(argv[i]);
    }

    std::unique_ptr<SpecularCoordinatesBrdf> outBrdf(DdrWriter::convert(*inBrdf));

    if (arranged) {
        outBrdf.reset(DdrWriter::arrange(*outBrdf, dataType));
    }

    // Save a DDR or DDT file.
    if (DdrWriter::write(outFileName, *outBrdf, comments)) {
        cout << "Converted: " << outFileName << endl;
        return 0;
    }
    else {
        cerr << "Failed to save: " << outFileName << endl;
        return 1;
    }
}
