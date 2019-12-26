// =================================================================== //
// Copyright (C) 2019 Kimura Ryo                                       //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#include <iostream>
#include <memory>

#include <libbsdf/Brdf/Analyzer.h>
#include <libbsdf/Common/SphericalCoordinateSystem.h>
#include <libbsdf/Reader/ReaderUtility.h>

#include <ArgumentParser.h>
#include <Utility.h>

using namespace lb;

/*
 * BRDF/BTDF analyzer
 */

const std::string APP_NAME("lbanalyze");
const std::string APP_VERSION("1.0.0");

const std::string ValueName         = "value";
const std::string ReflectanceName   = "reflectance";

// Paramters
float incomingPolarAngle = 0;
float incomingAzimuthalAngle = 0;
float outgoingPolarAngle = 0;
float outgoingAzimuthalAngle = 0;

void showHelp()
{
    using std::cout;
    using std::endl;

    cout << "Usage: lbanalyze [options ...] attribute_type file" << endl;
    cout << endl;
    cout << "lbanalyze analyzes a BRDF/BTDF and displays the result." << endl;
    cout << endl;
    cout << "Positional Arguments:" << endl;
    cout << "  attribute_type   Type of the analyzed attribute" << endl;
    cout << "  file             Name of the input BRDF/BTDF file." << endl;
    cout << "                   Valid formats:" << endl;
    cout << "                       Integra Diffuse Distribution (\".ddr, .ddt\")" << endl;
    cout << "                       LightTools/Zemax BSDF (\".bsdf\")" << endl;
    cout << "                       ASTM E1392-96(2002) (\".astm\")" << endl;
    cout << "                       MERL binary Files (\".binary\")" << endl;
    cout << endl;
    cout << "Options:" << endl;
    cout << "  -h, --help       show this help message and exit" << endl;
    cout << "  -v, --version    show program's version number and exit" << endl;
    cout << "  -l, --list       show acceptable attribute types and exit" << endl;
}

void showList()
{
    using std::cout;
    using std::endl;

    cout << "Acceptable attribute types:" << endl;
    cout << "  " << ValueName << endl;
    cout << "      Sampling of the BRDF/BTDF value at specified incoming and outgoing directions." << endl;
    cout << "      Valid options:" << endl;
    cout << "          -incomingPolarAngle      set the incoming polar angle in degrees (default: 0, range: [0, 90])" << endl;
    cout << "          -incomingAzimuthalAngle  set the incoming azimuthal angle in degrees (default: 0, range: [0, 360])" << endl;
    cout << "          -outgoingPolarAngle      set the outgoing polar angle in degrees (default: 0, range: [0, 90])" << endl;
    cout << "          -outgoingAzimuthalAngle  set the outgoing azimuthal angle in degrees (default: 0, range: [0, 360])" << endl;
    cout << "  " << ReflectanceName << endl;
    cout << "      Computation of the reflectance at a specified incoming direction." << endl;
    cout << "      Valid options:" << endl;
    cout << "          -incomingPolarAngle      set the incoming polar angle in degrees (default: 0, range: [0, 90])" << endl;
    cout << "          -incomingAzimuthalAngle  set the incoming azimuthal angle in degrees (default: 0, range: [0, 360])" << endl;
}

bool readOptions(ArgumentParser* ap)
{
    if (ap->read("-incomingPolarAngle", &incomingPolarAngle) == ArgumentParser::ERROR) {
        return false;
    }

    if (ap->read("-incomingAzimuthalAngle", &incomingAzimuthalAngle) == ArgumentParser::ERROR) {
        return false;
    }

    if (ap->read("-outgoingPolarAngle", &outgoingPolarAngle) == ArgumentParser::ERROR) {
        return false;
    }

    if (ap->read("-outgoingAzimuthalAngle", &outgoingAzimuthalAngle) == ArgumentParser::ERROR) {
        return false;
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

    if (ap.read("-l") || ap.read("--list")) {
        showList();
        return 0;
    }

    if (!readOptions(&ap)) return 1;

    if (!ap.validateNumTokens(2)) return 1;

    std::string attributeTypeName = ap.getTokens().at(0);
    std::string fileName          = ap.getTokens().at(1);

    // Load a BRDF/BTDF.
    FileType fileType;
    DataType dataType;
    std::shared_ptr<Brdf> brdf(reader_utility::read(fileName, &fileType, &dataType));

    if (!brdf) {
        std::cerr << "Failed to load: " << fileName << std::endl;
        return 1;
    }

    if (!brdf->validate()) {
        std::cerr << "Invalid attributes are found." << std::endl;
        return 1;
    }

    // Get an incoming direction.
    float inTheta = toRadian(incomingPolarAngle);
    float inPhi   = toRadian(incomingAzimuthalAngle);
    Vec3 inDir = SphericalCoordinateSystem::toXyz(inTheta, inPhi);

    // Get an outgoing direction.
    float outTheta = toRadian(outgoingPolarAngle);
    float outPhi   = toRadian(outgoingAzimuthalAngle);
    Vec3 outDir = SphericalCoordinateSystem::toXyz(outTheta, outPhi);

    if (attributeTypeName == ValueName) {
        // Sample and display the value of BRDF/BTDF.
        Spectrum value = brdf->getSpectrum(inDir, outDir);
        std::cout << "Value: " << value.format(LB_EIGEN_IO_FMT) << std::endl;
    }
    else if (attributeTypeName == ReflectanceName) {
        // Compute and display a reflectance array.
        Spectrum reflectance = computeReflectance(*brdf, inDir);
        std::cout << "Reflectance: " << reflectance.format(LB_EIGEN_IO_FMT) << std::endl;
    }

    return 0;
}
