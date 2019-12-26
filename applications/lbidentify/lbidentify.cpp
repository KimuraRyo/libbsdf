// =================================================================== //
// Copyright (C) 2019 Kimura Ryo                                       //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#include <iostream>
#include <memory>

#include <libbsdf/Brdf/HalfDifferenceCoordinatesBrdf.h>
#include <libbsdf/Brdf/SpecularCoordinatesBrdf.h>
#include <libbsdf/Brdf/SphericalCoordinatesBrdf.h>
#include <libbsdf/Reader/ReaderUtility.h>

#include <ArgumentParser.h>
#include <Utility.h>

using namespace lb;
using std::cout;
using std::cerr;
using std::endl;

/*
 * BRDF/BTDF identifier
 */

const std::string APP_NAME("lbidentify");
const std::string APP_VERSION("1.0.2");

void showHelp()
{
    cout << "Usage: lbidentify [options ...] file" << endl;
    cout << endl;
    cout << "lbidentify identifies a BRDF/BTDF file." << endl;
    cout << endl;
    cout << "Positional Arguments:" << endl;
    cout << "  file     Name of the input BRDF/BTDF file." << endl;
    cout << "           Valid formats:" << endl;
    cout << "               Integra Diffuse Distribution (\".ddr, .ddt\")" << endl;
    cout << "               LightTools/Zemax BSDF (\".bsdf\")" << endl;
    cout << "               ASTM E1392-96(2002) (\".astm\")" << endl;
    cout << "               MERL binary Files (\".binary\")" << endl;
    cout << endl;
    cout << "Options:" << endl;
    cout << "  -h, --help       show this help message and exit" << endl;
    cout << "  -v, --version    show program's version number and exit" << endl;
}

void showFileType(FileType fileType)
{
    switch (fileType) {
        case ASTM_FILE:
            cout << "File type: ASTM E1392-96(2002)" << endl;
            break;
        case INTEGRA_DDR_FILE:
            cout << "File type: Integra Diffuse Distribution Reflection" << endl;
            break;
        case INTEGRA_DDT_FILE:
            cout << "File type: Integra Diffuse Distribution Transparent" << endl;
            break;
        case lb::LIGHTTOOLS_FILE:
            cout << "File type: LightTools BSDF" << endl;
            break;
        case lb::MERL_BINARY_FILE:
            cout << "File type: MERL BRDF" << endl;
            break;
        case lb::ZEMAX_FILE:
            cout << "File type: Zemax BSDF" << endl;
            break;
        default:
            cerr << "Unknown file type: " << fileType << endl;
            break;
    }
}

void showDataType(DataType dataType)
{
    switch (dataType) {
        case lb::BRDF_DATA:
            cout << "Data type: BRDF" << endl;
            break;
        case lb::BTDF_DATA:
            cout << "Data type: BTDF" << endl;
            break;
        case lb::UNKNOWN_DATA:
            cout << "Data type is not distinguished between BRDF and BTDF." << endl;
            break;
        default:
            cerr << "Unsupported data type: " << dataType << endl;
            break;
    }
}

void showSourceType(SourceType sourceType)
{
    switch (sourceType) {
        case lb::MEASURED_SOURCE:
            cout << "Source type: Measured" << endl;
            break;
        case lb::EDITED_SOURCE:
            cout << "Source type: Edited" << endl;
            break;
        case lb::GENERATED_SOURCE:
            cout << "Source type: Generated" << endl;
            break;
        case lb::UNKNOWN_SOURCE:
            cout << "Source type: Unknown" << endl;
            break;
        default:
            cerr << "Unsupported data type: " << sourceType << endl;
            break;
    }
}

void showColorInfo(const SampleSet& samples)
{
    ColorModel cm = samples.getColorModel();
    switch (cm) {
        case lb::MONOCHROMATIC_MODEL:
            cout << "Color model: Monochromatic" << endl;
            break;
        case lb::RGB_MODEL:
            cout << "Color model: RGB" << endl;
            break;
        case lb::XYZ_MODEL:
            cout << "Color model: XYZ" << endl;
            break;
        case lb::SPECTRAL_MODEL:
            cout << "Color model: Spectral" << endl;
            cout << "Number of wavelengths: " << samples.getNumWavelengths() << endl;
            cout << "Wavelength (nm): " << samples.getWavelengths().format(LB_EIGEN_IO_FMT) << endl;
            break;
        default:
            cerr << "Unknown color model: " << cm << endl;
            break;
    }
}

void showAngleInfo(const Brdf& brdf)
{
    auto halfDiffBrdf = dynamic_cast<const HalfDifferenceCoordinatesBrdf*>(&brdf);
    auto specBrdf     = dynamic_cast<const SpecularCoordinatesBrdf*>(&brdf);
    auto spheBrdf     = dynamic_cast<const SphericalCoordinatesBrdf*>(&brdf);

    if (halfDiffBrdf) {
        cout << "Type of parameterization: Half difference" << endl;
    }
    else if (specBrdf) {
        cout << "Type of parameterization: Specular" << endl;
    }
    else if (spheBrdf) {
        cout << "Type of parameterization: Spherical" << endl;
    }
    else {
        cerr << "Failded to downcast lb::Brdf." << endl;
        return;
    }

    const SampleSet* ss = brdf.getSampleSet();

    using reader_utility::toLower;

    cout << "Number of " << toLower(brdf.getAngle0Name()) << "s: " << ss->getNumAngles0() << endl;
    cout << brdf.getAngle0Name() << ": " << toDegrees(ss->getAngles0()).format(LB_EIGEN_IO_FMT) << endl;

    cout << "Number of " << toLower(brdf.getAngle1Name()) << "s: " << ss->getNumAngles1() << endl;
    cout << brdf.getAngle1Name() << ": " << toDegrees(ss->getAngles1()).format(LB_EIGEN_IO_FMT) << endl;

    cout << "Number of " << toLower(brdf.getAngle2Name()) << "s: " << ss->getNumAngles2() << endl;
    cout << brdf.getAngle2Name() << ": " << toDegrees(ss->getAngles2()).format(LB_EIGEN_IO_FMT) << endl;

    cout << "Number of " << toLower(brdf.getAngle3Name()) << "s: " << ss->getNumAngles3() << endl;
    cout << brdf.getAngle3Name() << ": " << toDegrees(ss->getAngles3()).format(LB_EIGEN_IO_FMT) << endl;

    if (specBrdf && specBrdf->getNumSpecularOffsets() != 0) {
        cout << "Number of specular offsets: " << specBrdf->getNumSpecularOffsets() << endl;
        cout << "Specular offsets: " << toDegrees(specBrdf->getSpecularOffsets()).format(LB_EIGEN_IO_FMT) << endl;
    }
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

    if (!ap.validateNumTokens(1)) return 1;

    std::string fileName = ap.getTokens().at(0);

    // Load a BRDF/BTDF.
    FileType fileType;
    DataType dataType;
    std::shared_ptr<Brdf> brdf(reader_utility::read(fileName, &fileType, &dataType));

    if (brdf) {
        cout << "File name: " << fileName << endl;
    }
    else {
        cerr << "Failed to load: " << fileName << endl;
        return 1;
    }

    if (!brdf->validate()) {
        cout << "Invalid attributes are found." << endl;
    }

    // Display information.
    showFileType(fileType);
    showDataType(dataType);
    showSourceType(brdf->getSourceType());
    showAngleInfo(*brdf);
    showColorInfo(*brdf->getSampleSet());

    return 0;
}
