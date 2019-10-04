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
#include <libbsdf/Brdf/SpecularCoordinatesBrdf.h>

#include <libbsdf/Common/Log.h>

#include <libbsdf/Reader/DdrReader.h>
#include <libbsdf/Reader/ReaderUtility.h>

#include <libbsdf/Writer/DdrWriter.h>

#include <ArgumentParser.h>
#include <Utility.h>

using namespace lb;

/*
 * BRDF/BTDF comparer
 */

const std::string APP_NAME("lbcompare");
const std::string APP_VERSION("1.0.0");

void showHelp()
{
    using std::cout;
    using std::endl;

    cout << "Usage: lbcompare [options ...] in_file1 in_file2 out_file" << endl;
    cout << endl;
    cout << "lbcompare compares two BRDFs/BTDFs and saves the difference between them." << endl;
    cout << "Two input files must have the same type of color." << endl;
    cout << endl;
    cout << "Positional Arguments:" << endl;
    cout << "  in_file1, in_file2   Names of input BRDF/BTDF files." << endl;
    cout << "                       Valid formats:" << endl;
    cout << "                           Integra Diffuse Distribution (\".ddr, .ddt\")" << endl;
    cout << "  out_file             Name of an output BRDF/BTDF file. Attributes of angles are same as in_file1." << endl;
    cout << "                       Valid formats:" << endl;
    cout << "                           Integra Diffuse Distribution (\".ddr, .ddt\")" << endl;
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

    std::string inFileName1 = ap.getTokens().at(0);
    std::string inFileName2 = ap.getTokens().at(1);
    std::string outFileName = ap.getTokens().at(2);

    // Validate file types.
    FileType inFileType = reader_utility::classifyFile(inFileName1);
    if (inFileType != reader_utility::classifyFile(inFileName2)) {
        std::cerr << "Input file types do not match." << std::endl;
        return 1;
    }
    else if (inFileType != INTEGRA_DDR_FILE &&
             inFileType != INTEGRA_DDT_FILE) {
        std::cerr << "Unsupported file type: " << inFileType << std::endl;
        return 1;
    }

    // Load BRDF/BTDF files.
    std::unique_ptr<SpecularCoordinatesBrdf> inBrdf1(DdrReader::read(inFileName1));
    std::unique_ptr<SpecularCoordinatesBrdf> inBrdf2(DdrReader::read(inFileName2));

    if (!inBrdf1) {
        std::cerr << "Failed to load: " << inFileName1 << std::endl;
        return 1;
    }
    else if (!inBrdf2) {
        std::cerr << "Failed to load: " << inFileName2 << std::endl;
        return 1;
    }

    // Validate color models.
    if (!hasSameColor(*inBrdf1->getSampleSet(), *inBrdf2->getSampleSet())) {
        std::cerr << "Color models or wavelengths do not match." << std::endl;
        return 1;
    }

    // Make a BRDF with the absolute difference between two inputs.
    std::unique_ptr<SpecularCoordinatesBrdf> outBrdf(inBrdf1->clone());
    subtract(*inBrdf1, *inBrdf2, outBrdf.get());
    for (Spectrum& sp : outBrdf->getSampleSet()->getSpectra()) {
        sp = sp.cwiseAbs();
    }

    // Save a DDR/DDT file.
    std::string comments = app_utility::createComments(argc, argv, APP_NAME, APP_VERSION);
    if (DdrWriter::write(outFileName, *outBrdf, comments)) {
        std::cout << "Saved: " << outFileName << std::endl;
    }

    return 0;
}
