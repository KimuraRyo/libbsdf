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

#include <libbsdf/Reader/DdrReader.h>
#include <libbsdf/Reader/ReaderUtility.h>

#include <libbsdf/Writer/DdrWriter.h>

#include <ArgumentParser.h>
#include <Utility.h>

using namespace lb;

/*
 * BRDF/BTDF averager
 */

const std::string APP_NAME("lbaverage");
const std::string APP_VERSION("1.0.0");

void showHelp()
{
    using std::cout;
    using std::endl;

    cout << "Usage: lbaverage [options ...] in_file1 in_file2 ... out_file" << endl;
    cout << endl;
    cout << "lbaverage averages two or more BRDFs/BTDFs." << endl;
    cout << "Input files must have the same type of color." << endl;
    cout << endl;
    cout << "Positional Arguments:" << endl;
    cout << "  in_file1, in_file2...    Names of input BRDF/BTDF files." << endl;
    cout << "                           Valid formats:" << endl;
    cout << "                               Integra Diffuse Distribution (\".ddr, .ddt\")" << endl;
    cout << "  out_file                 Name of an output BRDF/BTDF file. Attributes of angles are same as in_file1." << endl;
    cout << "                           Valid formats:" << endl;
    cout << "                               Integra Diffuse Distribution (\".ddr, .ddt\")" << endl;
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

    const std::vector<std::string>& fileNames = ap.getTokens();
    int numInputFiles = static_cast<int>(fileNames.size()) - 1;

    if (numInputFiles < 2) {
        std::cerr << "Two or more input files are required." << std::endl;
        return 1;
    }

    const std::string& inFileName1 = fileNames.at(0);

    // Set and validate the file type of the first input file.
    FileType inFileType = reader_utility::classifyFile(inFileName1);
    if (inFileType != INTEGRA_DDR_FILE &&
        inFileType != INTEGRA_DDT_FILE) {
        std::cerr << "Unsupported file type: " << inFileType << std::endl;
        return 1;
    }

    // Load the first input BRDF/BTDF file.
    std::unique_ptr<SpecularCoordinatesBrdf> inBrdf1(DdrReader::read(inFileName1));
    if (!inBrdf1) {
        std::cerr << "Failed to load: " << inFileName1 << std::endl;
        return 1;
    }

    // Create an output BRDF with the same angle attributes as in_file1;
    std::unique_ptr<SpecularCoordinatesBrdf> outBrdf(inBrdf1->clone());

    Spectrum::Scalar avgCoeff = Spectrum::Scalar(1) / numInputFiles;

    // Set a scaled BRDF to compute the average of spectra.
    multiplySpectra(outBrdf->getSampleSet(), avgCoeff);

    // Read input files and compute averaged spectra.
    for (size_t i = 1; i < numInputFiles; ++i) {
        const std::string& fileName = fileNames.at(i);

        // Validate a file type.
        if (inFileType != reader_utility::classifyFile(fileName)) {
            std::cerr << "Input file types do not match: " << fileName << std::endl;
            return 1;
        }

        // Load an input BRDF/BTDF file.
        std::unique_ptr<SpecularCoordinatesBrdf> inBrdf(DdrReader::read(fileName));
        if (!inBrdf) {
            std::cerr << "Failed to load: " << fileName << std::endl;
            return 1;
        }

        // Copy previously processed BRDF.
        std::unique_ptr<SpecularCoordinatesBrdf> prevOutBrdf(outBrdf->clone());

        // Add a scaled BRDF to compute the average of spectra.
        auto add = [avgCoeff](const Spectrum& sp1, const Spectrum& sp2) { return sp1 + sp2 * avgCoeff; };
        if (!compute(*prevOutBrdf, *inBrdf, outBrdf.get(), add)) {
            std::cerr << "Failed to process: " << fileName << std::endl;
            return 1;
        }
    }

    // Save a DDR/DDT file.
    std::string outFileName = fileNames.back();
    std::string comments = app_utility::createComments(argc, argv, APP_NAME, APP_VERSION);
    if (DdrWriter::write(outFileName, *outBrdf, comments)) {
        std::cout << "Saved: " << outFileName << std::endl;
    }

    return 0;
}
