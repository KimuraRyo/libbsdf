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

#include <libbsdf/Writer/DdrWriter.h>

#include <ArgumentParser.h>
#include <Utility.h>

using namespace lb;

/*
 * BRDF inserter
 */

int main(int argc, char** argv)
{
    Log::setNotificationLevel(Log::Level::WARN_MSG);

    ArgumentParser ap(argc, argv);

    using std::cout;
    using std::cerr;
    using std::endl;

    if (ap.read("-h") || ap.read("--help") || ap.getTokens().empty()) {
        cout << "Usage: lbinsert [options ...] partial_file angle base_file" << endl;
        cout << endl;
        cout << "lbinsert inserts a BRDF with an incoming azimuthal angle into a base BRDF." << endl;
        cout << "Two files must have the same type of color." << endl;
        cout << endl;
        cout << "Positional Arguments:" << endl;
        cout << "  partial_file Name of a BRDF file inserted into a base BRDF." << endl;
        cout << "               This BRDF must have an incoming azimuthal angle." << endl;
        cout << "               \".ddr\" or \".astm\" is acceptable as a suffix." << endl;
        cout << "  angle        Incoming azimuthal angle in degrees (range: [0, 360])" << endl;
        cout << "  base_file    Name of a base BRDF file." << endl;
        cout << "               \".ddr\" is acceptable as a suffix." << endl;
        cout << endl;
        cout << "Options:" << endl;
        cout << "  -h, --help       show this help message and exit" << endl;
        cout << "  -v, --version    show program's version number and exit" << endl;
        return 0;
    }

    const std::string version("1.0.0");

    if (ap.read("-v") || ap.read("--version")) {
        cout << "Version: lbinsert " << version << " (libbsdf-" << getVersion() << ")" << endl;
        return 0;
    }

    if (!ap.validate(3)) return 1;

    std::string partialFileName = ap.getTokens().at(0);
    std::string angleStr        = ap.getTokens().at(1);
    std::string baseFileName    = ap.getTokens().at(2);

    FileType partialFileType = reader_utility::classifyFile(partialFileName);
    FileType baseFileType    = reader_utility::classifyFile(baseFileName);

    std::unique_ptr<SpecularCoordinatesBrdf> baseBrdf;

    // Load a base BRDF file.
    switch (baseFileType) {
        case INTEGRA_DDR_FILE:
            baseBrdf.reset(DdrReader::read(baseFileName));
            break;
        default:
            cerr << "Invalid file type: " << baseFileName << endl;
            return 1;
    }

    if (!baseBrdf) {
        cerr << "Failed to load: " << baseFileName << endl;
        return 1;
    }

    std::unique_ptr<SpecularCoordinatesBrdf> partialBrdf;

    // Load a partial BRDF file.
    switch (partialFileType) {
        case ASTM_FILE: {
            std::unique_ptr<SphericalCoordinatesBrdf> brdf(AstmReader::read(partialFileName));
            if (!brdf) {
                cerr << "Failed to load: " << partialFileName << endl;
                return 1;
            }

            partialBrdf.reset(new SpecularCoordinatesBrdf(*brdf, baseBrdf->getNumSpecTheta(), baseBrdf->getNumSpecPhi()));
            break;
        }
        case INTEGRA_DDR_FILE:
            partialBrdf.reset(DdrReader::read(partialFileName));
            break;
        default:
            cerr << "Invalid file type: " << partialFileName << endl;
            return 1;
    }

    if (!partialBrdf) {
        cerr << "Failed to load: " << partialFileName << endl;
        return 1;
    }

    SampleSet* baseSs    = baseBrdf->getSampleSet();
    SampleSet* partialSs = partialBrdf->getSampleSet();

    // If two monochromatic data sets have different color modes, partial data is adjusted.
    if (!hasSameColor(*baseSs, *partialSs) &&
        partialSs->getColorModel() == SPECTRAL_MODEL &&
        partialSs->getNumWavelengths() == 1) {
        partialSs->setColorModel(MONOCHROMATIC_MODEL);
        partialSs->setWavelength(0, 0.0f);
    }

    // Read an incoming azimuthal angle.
    char* end;
    double inPhiDegree = std::strtod(angleStr.c_str(), &end);
    if (*end != '\0') {
        cerr << "Invalid value: " << angleStr << endl;
        return 1;
    }

    inPhiDegree = utility::clampParameter("angle", inPhiDegree, 0.0, 360.0);

    float inPhi = static_cast<float>(lb::toRadian(inPhiDegree));
    std::unique_ptr<SpecularCoordinatesBrdf> brdf(insertBrdfAlongInPhi(*baseBrdf, *partialBrdf, inPhi));

    std::string comments("Software: lbinsert-" + version);
    comments += "\n;; Arguments:";
    for (int i = 1; i < argc; ++i) {
        comments += " " + std::string(argv[i]);
    }

    DdrWriter::write(baseFileName, *brdf, comments);

    return 0;
}
