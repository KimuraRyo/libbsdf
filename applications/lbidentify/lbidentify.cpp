// =================================================================== //
// Copyright (C) 2019-2020 Kimura Ryo                                  //
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
const std::string APP_VERSION("1.0.3");

void showHelp()
{
    cout << "Usage: lbidentify [options ...] file" << endl;
    cout << endl;
    cout << "lbidentify identifies a BRDF/BTDF file." << endl;
    cout << endl;
    cout << "Positional Arguments:" << endl;
    cout << "  file     Name of the input BRDF/BTDF file." << endl;
    cout << "           Valid formats:" << endl;
    cout << "               Surface Scattering Distribution Data (\".ssdd\")" << endl;
    cout << "               Integra Diffuse Distribution (\".ddr, .ddt\")" << endl;
    cout << "               LightTools/Zemax (\".bsdf\")" << endl;
    cout << "               ASTM E1392-96(2002) (\".astm\")" << endl;
    cout << "               MERL binary (\".binary\")" << endl;
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
        case LIGHTTOOLS_FILE:
            cout << "File type: LightTools BSDF" << endl;
            break;
        case MERL_BINARY_FILE:
            cout << "File type: MERL BRDF" << endl;
            break;
        case SSDD_FILE:
            cout << "File type: Surface Scattering Distribution Data" << endl;
            break;
        case ZEMAX_FILE:
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
        case BRDF_DATA:
            cout << "Data type: BRDF" << endl;
            break;
        case BTDF_DATA:
            cout << "Data type: BTDF" << endl;
            break;
        case SPECULAR_REFLECTANCE_DATA:
            cout << "Data type: Specular reflectance" << endl;
            break;
        case SPECULAR_TRANSMITTANCE_DATA:
            cout << "Data type: Specular transmittance" << endl;
            break;
        case UNKNOWN_DATA:
            cout << "Data type: unknown" << endl;
            break;
        default:
            cerr << "Unsupported data type: " << dataType << endl;
            break;
    }
}

void showSourceType(SourceType sourceType)
{
    switch (sourceType) {
        case MEASURED_SOURCE:
            cout << "Source type: Measured" << endl;
            break;
        case EDITED_SOURCE:
            cout << "Source type: Edited" << endl;
            break;
        case GENERATED_SOURCE:
            cout << "Source type: Generated" << endl;
            break;
        case UNKNOWN_SOURCE:
            cout << "Source type: Unknown" << endl;
            break;
        default:
            cerr << "Unsupported data type: " << sourceType << endl;
            break;
    }
}

void showAngleInfo(const Brdf& brdf)
{
    auto halfDiffBrdf = dynamic_cast<const HalfDifferenceCoordinatesBrdf*>(&brdf);
    auto specBrdf     = dynamic_cast<const SpecularCoordinatesBrdf*>(&brdf);
    auto spheBrdf     = dynamic_cast<const SphericalCoordinatesBrdf*>(&brdf);

    if (halfDiffBrdf) {
        cout << "Type of parameterization: Half difference coordinate system" << endl;
    }
    else if (specBrdf) {
        cout << "Type of parameterization: Specular coordinate system" << endl;
    }
    else if (spheBrdf) {
        cout << "Type of parameterization: Spherical coordinate system" << endl;
    }
    else {
        cerr << "Failded to downcast lb::Brdf." << endl;
        return;
    }

    const SampleSet* ss = brdf.getSampleSet();

    using reader_utility::toLower;

    cout << "  Name: " << brdf.getAngle0Name() << endl;
    cout << "    Size: " << ss->getNumAngles0() << endl;
    cout << "    Degree: " << toDegrees(ss->getAngles0()).format(LB_EIGEN_IO_FMT) << endl;

    cout << "  Name: " << brdf.getAngle1Name() << endl;
    cout << "    Size: " << ss->getNumAngles1() << endl;
    cout << "    Degree: " << toDegrees(ss->getAngles1()).format(LB_EIGEN_IO_FMT) << endl;

    cout << "  Name: " << brdf.getAngle2Name() << endl;
    cout << "    Size: " << ss->getNumAngles2() << endl;
    cout << "    Degree: " << toDegrees(ss->getAngles2()).format(LB_EIGEN_IO_FMT) << endl;

    cout << "  Name: " << brdf.getAngle3Name() << endl;
    cout << "    Size: " << ss->getNumAngles3() << endl;
    cout << "    Degree: " << toDegrees(ss->getAngles3()).format(LB_EIGEN_IO_FMT) << endl;

    if (specBrdf && specBrdf->getNumSpecularOffsets() != 0) {
        cout << "  Name: Specular offset" << endl;
        cout << "    Size: " << specBrdf->getNumSpecularOffsets() << endl;
        cout << "    Degree: " << toDegrees(specBrdf->getSpecularOffsets()).format(LB_EIGEN_IO_FMT) << endl;
    }
}

void showAngleInfo(const SampleSet2D& ss2)
{
    cout << "Type of parameterization: Spherical coordinate system" << endl;
    cout << "  Name: Incoming polar angle" << endl;
    cout << "    Size: " << ss2.getNumTheta() << endl;
    cout << "    Degree: " << toDegrees(ss2.getThetaArray()).format(LB_EIGEN_IO_FMT) << endl;

    cout << "  Name: Incoming azimuthal angle" << endl;
    cout << "    Size: " << ss2.getNumPhi() << endl;
    cout << "    Degree: " << toDegrees(ss2.getPhiArray()).format(LB_EIGEN_IO_FMT) << endl;
}

// Show color information of BRDF, BTDF, specular reflectance, or specular transmittance.
template <typename T>
void showColorInfo(const T& data)
{
    ColorModel cm = data.getColorModel();
    switch (cm) {
        case MONOCHROMATIC_MODEL:
            cout << "Color model: Monochrome" << endl;
            break;
        case RGB_MODEL:
            cout << "Color model: RGB" << endl;
            break;
        case XYZ_MODEL:
            cout << "Color model: XYZ" << endl;
            break;
        case SPECTRAL_MODEL:
            cout << "Color model: Spectrum" << endl;
            cout << "  Size: " << data.getNumWavelengths() << endl;
            cout << "  Wavelength (nm): " << data.getWavelengths().format(LB_EIGEN_IO_FMT) << endl;
            break;
        default:
            cerr << "Unknown color model: " << cm << endl;
            break;
    }
}

void showDataInfo(const Brdf& brdf, DataType dataType)
{
    if (!brdf.validate()) {
        cout << "Invalid attributes are found." << endl;
        return;
    }

    cout << endl;
    showDataType(dataType);
    showSourceType(brdf.getSourceType());
    showAngleInfo(brdf);
    showColorInfo(*brdf.getSampleSet());
}

void showDataInfo(const SampleSet2D& ss2, DataType dataType)
{
    if (!ss2.validate()) {
        cout << "Invalid attributes are found." << endl;
        return;
    }

    cout << endl;
    showDataType(dataType);
    showSourceType(ss2.getSourceType());
    showAngleInfo(ss2);
    showColorInfo(ss2);
}

void showBrdfInfo(const std::string& fileName)
{
    FileType fileType;
    DataType dataType;
    std::shared_ptr<Brdf> brdf = reader_utility::readBrdf(fileName, &fileType, &dataType);

    if (brdf) {
        cout << "File name: " << fileName << endl;
    }
    else {
        cerr << "Failed to load: " << fileName << endl;
        return;
    }

    showDataInfo(*brdf, dataType);
}

void showMaterialInfo(const std::string& fileName)
{
    FileType fileType;
    std::shared_ptr<Material> material = reader_utility::readMaterial(fileName, &fileType);

    if (material) {
        cout << "File name: " << fileName << endl;
    }
    else {
        cerr << "Failed to load: " << fileName << endl;
        return;
    }

    if (std::shared_ptr<Bsdf> bsdf = material->getBsdf()) {
        if (std::shared_ptr<Brdf> brdf = bsdf->getBrdf()) {
            showDataInfo(*brdf, BRDF_DATA);
        }

        if (std::shared_ptr<Btdf> btdf = bsdf->getBtdf()) {
            showDataInfo(*btdf->getBrdf(), BTDF_DATA);
        }
    }

    if (std::shared_ptr<SampleSet2D> sr = material->getSpecularReflectances()) {
        showDataInfo(*sr, SPECULAR_REFLECTANCE_DATA);
    }

    if (std::shared_ptr<SampleSet2D> st = material->getSpecularTransmittances()) {
        showDataInfo(*st, SPECULAR_TRANSMITTANCE_DATA);
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

    FileType fileType = reader_utility::classifyFile(fileName);
    showFileType(fileType);

    // Load data and show information.
    switch (fileType) {
        case LIGHTTOOLS_FILE:
        case SSDD_FILE:
            showMaterialInfo(fileName);
            break;
        case INTEGRA_SDR_FILE:
        case INTEGRA_SDT_FILE:
            cout << "Unsupported file type." << endl;
            return 0;
        case UNKNOWN_FILE:
            cerr << "Unknown file type." << endl;
            return 1;
        default:
            showBrdfInfo(fileName);
            break;
    }

    return 0;
}
