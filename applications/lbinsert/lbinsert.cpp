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

#include <libbsdf/Common/Version.h>

#include <libbsdf/Reader/AstmReader.h>
#include <libbsdf/Reader/DdrReader.h>
#include <libbsdf/Reader/ReaderUtility.h>

#include <libbsdf/Writer/DdrWriter.h>

#include <ArgumentParser.h>

using namespace lb;

/*
 * BRDF inserter
 */

int main(int argc, char** argv)
{
    ArgumentParser ap(argc, argv);

    if (ap.read("-h") || ap.read("--help") || ap.getTokens().empty()) {
        std::cout << "Usage: lbinsert [options ...] partial_file angle base_file" << std::endl;
        std::cout << std::endl;
        std::cout << "lbinsert inserts a BRDF with an incoming azimuthal angle into a base BRDF." << std::endl;
        std::cout << "Two files must have the same type of color." << std::endl;
        std::cout << std::endl;
        std::cout << "Positional Arguments:" << std::endl;
        std::cout << "  partial_file    Name of a BRDF file inserted into a base BRDF." << std::endl;
        std::cout << "                  This BRDF must have an incoming azimuthal angle." << std::endl;
        std::cout << "                  \".ddr\" or \".astm\" is acceptable as a suffix." << std::endl;
        std::cout << "  angle           Incoming azimuthal angle in degrees" << std::endl;
        std::cout << "  base_file       Name of a base BRDF file." << std::endl;
        std::cout << "                  \".ddr\" is acceptable as a suffix." << std::endl;
        std::cout << std::endl;
        std::cout << "Options:" << std::endl;
        std::cout << "  -h, --help      show this help message and exit" << std::endl;
        std::cout << "  -v, --version   show program's version number and exit" << std::endl;
        return 0;
    }

    const std::string version("1.0.0");

    if (ap.read("-v") || ap.read("--version")) {
        std::cout << "Version: lbinsert " << version << " (libbsdf-" << getVersion() << ")" << std::endl;
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
            std::cerr << "Invalid file type: " << baseFileName << std::endl;
            return 1;
    }

    if (!baseBrdf) {
        std::cerr << "Failed to load: " << baseFileName << std::endl;
        return 1;
    }

    std::unique_ptr<SpecularCoordinatesBrdf> partialBrdf;

    // Load a partial BRDF file.
    switch (partialFileType) {
        case ASTM_FILE: {
            std::unique_ptr<SphericalCoordinatesBrdf> brdf(AstmReader::read(partialFileName));
            if (!brdf) {
                std::cerr << "Failed to load: " << partialFileName << std::endl;
                return 1;
            }

            partialBrdf.reset(new SpecularCoordinatesBrdf(*brdf, baseBrdf->getNumSpecTheta(), baseBrdf->getNumSpecPhi()));
            break;
        }
        case INTEGRA_DDR_FILE:
            partialBrdf.reset(DdrReader::read(partialFileName));
            break;
        default:
            std::cerr << "Invalid file type: " << partialFileName << std::endl;
            return 1;
    }

    if (!partialBrdf) {
        std::cerr << "Failed to load: " << partialFileName << std::endl;
        return 1;
    }

    // Read an incoming azimuthal angle.
    char* end;
    double inPhiDegree = std::strtod(angleStr.c_str(), &end);
    if (*end != '\0') {
        std::cerr << "Invalid value: " << angleStr << std::endl;
        return 1;
    }

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
