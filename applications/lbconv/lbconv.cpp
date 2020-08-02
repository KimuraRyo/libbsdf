// =================================================================== //
// Copyright (C) 2019-2020 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#include <iostream>
#include <memory>

#include <libbsdf/Reader/AstmReader.h>
#include <libbsdf/Reader/DdrReader.h>
#include <libbsdf/Reader/SsddReader.h>
#include <libbsdf/Reader/ZemaxBsdfReader.h>

#include <libbsdf/Writer/DdrWriter.h>
#include <libbsdf/Writer/SsddWriter.h>

#include <ArgumentParser.h>
#include <Utility.h>

using namespace lb;

/*
 * BRDF/BTDF converter
 */

const std::string APP_NAME("lbconv");
const std::string APP_VERSION("1.0.2");

// Parameters
DataType dataType = BRDF_DATA;
bool arranged = false;

void showHelp()
{
    using std::cout;
    using std::endl;

    cout << "Usage: lbconv [options ...] in_file out_file" << endl;
    cout << endl;
    cout << "lbconv converts a BRDF/BTDF file to another format." << endl;
    cout << endl;
    cout << "Positional Arguments:" << endl;
    cout << "  in_file      Name of an input BRDF/BTDF file." << endl;
    cout << "               Valid formats:" << endl;
    cout << "                   Surface Scattering Distribution Data (\".ssdd\")" << endl;
    cout << "                   Integra Diffuse Distribution (\".ddr, .ddt\")" << endl;
    cout << "                   Zemax BSDF (\".bsdf\")" << endl;
    cout << "                   ASTM E1392-96(2002) (\".astm\")" << endl;
    cout << "  out_file     Name of an output BRDF/BTDF file." << endl;
    cout << "               Valid formats:" << endl;
    cout << "                   Surface Scattering Distribution Data (\".ssdd\")" << endl;
    cout << "                   Integra Diffuse Distribution (\".ddr, .ddt\")" << endl;
    cout << "                   If an appropriate suffix is not obtained, \".ddr\" or \".ddt\" is appended." << endl;
    cout << endl;
    cout << "Options:" << endl;
    cout << "  -h, --help       show this help message and exit" << endl;
    cout << "  -v, --version    show program's version number and exit" << endl;
    cout << "  -scatterType     set either BRDF or BTDF for the input ASTM file (default: BRDF)" << endl;
    cout << "  -arrangement     arrange BRDF/BTDF with extrapolation and conservation of energy" << endl;
}

bool readOptions(ArgumentParser* ap)
{
    std::string dataTypeStr;
    if (ap->read("-scatterType", &dataTypeStr)) {
        if (dataTypeStr == "BRDF") {
            dataType = BRDF_DATA;
        }
        else if (dataTypeStr == "BTDF") {
            dataType = BTDF_DATA;
        }
        else {
            std::cerr << "Invalid scatter type: " << dataTypeStr << std::endl;
            return false;
        }
    }

    if (ap->read("-arrangement")) {
        arranged = true;
    }

    return true;
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

    if (!readOptions(&ap)) return 1;

    if (!ap.validateNumTokens(2)) return 1;

    std::string inFileName  = ap.getTokens().at(0);
    std::string outFileName = ap.getTokens().at(1);

    // Read a BRDF/BTDF file.
    FileType inFileType = reader_utility::classifyFile(inFileName);
    std::shared_ptr<Brdf> brdf;
    std::shared_ptr<Material> material;
    switch (inFileType) {
        case SSDD_FILE:
            material = reader_utility::readMaterial(inFileName, &inFileType);
            break;
        case ASTM_FILE:
        case INTEGRA_DDR_FILE:
        case INTEGRA_DDT_FILE:
        case ZEMAX_FILE:
            brdf = reader_utility::readBrdf(inFileName, &inFileType, &dataType);
            break;
        default:
            std::cerr << "Unsupported file type: " << inFileType << std::endl;
            return 1;
    }

    if (!brdf && !material) {
        std::cerr << "Failed to load: " << inFileName << std::endl;
        return 1;
    }

    if (material) {
        bool empty = (!material->getBsdf() || material->getBsdf()->isEmpty());
        if (empty) {
            std::cerr << "BRDF/BTDF is not found." << std::endl;
            return 1;
        }
    }

    std::string comments = app_utility::createComments(argc, argv, APP_NAME, APP_VERSION);

    if (reader_utility::hasSuffix(outFileName, ".ssdd")) {
        if (brdf) {
            std::shared_ptr<Bsdf> bsdf;
            if (dataType == BRDF_DATA) {
                bsdf = std::make_shared<Bsdf>(brdf, nullptr);
            }
            else if (dataType == BTDF_DATA) {
                std::shared_ptr<Btdf> btdf = std::make_shared<Btdf>(brdf);
                bsdf = std::make_shared<Bsdf>(nullptr, btdf);
            }
            else {
                std::cerr << "Invalid data type: " << dataType << std::endl;
                return 1;
            }

            material->setBsdf(bsdf);
        }

        // Write an SSDD file.
        SsddWriter::write(outFileName, *material, SsddWriter::DataFormat::ASCII_DATA, comments);
    }
    else {
        bool ddrSuffixFound = reader_utility::hasSuffix(outFileName, ".ddr");
        bool ddtSuffixFound = reader_utility::hasSuffix(outFileName, ".ddt");

        if (material) {
            std::shared_ptr<Bsdf> bsdf = material->getBsdf();

            if (bsdf->getBrdf() && ddrSuffixFound) {
                brdf = bsdf->getBrdf();
                dataType = BRDF_DATA;
            }
            else if (bsdf->getBtdf() && ddtSuffixFound) {
                brdf = bsdf->getBtdf()->getBrdf();
                dataType = BTDF_DATA;
            }
        }

        // Convert the BRDF/BTDF.
        std::unique_ptr<SpecularCoordinatesBrdf> outBrdf(DdrWriter::convert(*brdf));
        if (arranged) {
            outBrdf.reset(DdrWriter::arrange(*outBrdf, dataType));
        }

        // Fix the output filename.
        if (dataType == BRDF_DATA && !ddrSuffixFound) {
            outFileName += ".ddr";
        }
        else if (dataType == BTDF_DATA && !ddtSuffixFound) {
            outFileName += ".ddt";
        }

        // Write a DDR/DDT file.
        if (DdrWriter::write(outFileName, *outBrdf, comments)) {
            std::cout << "Saved: " << outFileName << std::endl;
        }
    }

    return 0;
}
