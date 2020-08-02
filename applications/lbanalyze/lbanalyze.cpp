// =================================================================== //
// Copyright (C) 2019-2020 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#include <iostream>
#include <memory>

#include <libbsdf/Brdf/Analyzer.h>
#include <libbsdf/Common/SpectrumUtility.h>
#include <libbsdf/Common/SphericalCoordinateSystem.h>
#include <libbsdf/Reader/ReaderUtility.h>

#include <ArgumentParser.h>
#include <Utility.h>

using namespace lb;

/*
 * BRDF/BTDF analyzer
 */

const std::string APP_NAME("lbanalyze");
const std::string APP_VERSION("1.0.2");

const std::string ValueName         = "value";
const std::string ReflectanceName   = "reflectance";

// Parameters
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
    cout << "                       Surface Scattering Distribution Data (\".ssdd\")" << endl;
    cout << "                       Integra Diffuse Distribution (\".ddr, .ddt\")" << endl;
    cout << "                       LightTools/Zemax (\".bsdf\")" << endl;
    cout << "                       ASTM E1392-96(2002) (\".astm\")" << endl;
    cout << "                       MERL binary (\".binary\")" << endl;
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
    cout << "          -outgoingPolarAngle      set the outgoing polar angle in degrees (default: 0, range: [0, 180])" << endl;
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

void showReflectance(const Vec3& xyz, const Vec3& rgb)
{
    Vec3 lab = xyzToLab(xyz);

    std::cout << "Linear sRGB: " << rgb.format(LB_EIGEN_IO_FMT) << std::endl;
    std::cout << "CIE XYZ: "     << xyz.format(LB_EIGEN_IO_FMT) << std::endl;
    std::cout << "CIE LAB: "     << lab.format(LB_EIGEN_IO_FMT) << std::endl;

    std::string hue;
    float value;
    int chroma;
    findMunsellProperties(xyz, &hue, &value, &chroma);
    std::cout << "Munsell: " << hue;
    if (chroma == 0.0f) {
        std::cout << value << std::endl;
    }
    else {
        std::cout << " " << value << "/" << chroma << std::endl;
    }
}

void showReflectance(const Spectrum& reflectance, ColorModel cm, const Arrayf& wavelengths, Vec3* xyz)
{
    switch (cm) {
        case MONOCHROMATIC_MODEL:
            std::cout << "Value: " << reflectance.format(LB_EIGEN_IO_FMT) << std::endl;
            break;
        case SPECTRAL_MODEL:
            std::cout << "Spectrum: " << reflectance.format(LB_EIGEN_IO_FMT) << std::endl;
            std::cout << "Wavelength (nm): " << wavelengths.format(LB_EIGEN_IO_FMT) << std::endl;
            break;
        default:
            break;
    }

    *xyz = SpectrumUtility::spectrumToXyz(reflectance, cm, wavelengths);

    Vec3 rgb = (cm == RGB_MODEL) ? toVec3(reflectance) : xyzToSrgb(*xyz).cwiseMax(0.0);
    showReflectance(*xyz, rgb);
}

void showReflectance(const Brdf* brdf, const Vec3& inDir, Vec3* xyz)
{
    const SampleSet* ss = brdf->getSampleSet();

    Spectrum reflectance = computeReflectance(*brdf, inDir);
    ColorModel cm = ss->getColorModel();
    const Arrayf& wavelengths = ss->getWavelengths();
    showReflectance(reflectance, cm, wavelengths, xyz);
}

void showReflectance(const SampleSet2D* specularReflectance, const Vec3& inDir, Vec3* xyz)
{
    Spectrum reflectance = specularReflectance->getSpectrum(inDir);
    ColorModel cm = specularReflectance->getColorModel();
    const Arrayf& wavelengths = specularReflectance->getWavelengths();
    showReflectance(reflectance, cm, wavelengths, xyz);
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

    // Load the material containing a BSDF.
    FileType fileType;
    std::shared_ptr<Material> material = reader_utility::readMaterial(fileName, &fileType);
    if (!material) {
        std::cerr << "Failed to load: " << fileName << std::endl;
        return 1;
    }

    std::shared_ptr<const Bsdf> bsdf = material->getBsdf();
    std::shared_ptr<const SampleSet2D> sr = material->getSpecularReflectances();
    std::shared_ptr<const SampleSet2D> st = material->getSpecularTransmittances();
    if ((bsdf && !bsdf->validate()) ||
        (sr   && !sr->validate()) ||
        (st   && !st->validate())) {
        std::cerr << "Invalid attributes are found." << std::endl;
        return 1;
    }

    // Get the incoming direction.
    float inTheta = toRadian(incomingPolarAngle);
    float inPhi   = toRadian(incomingAzimuthalAngle);
    Vec3 inDir = SphericalCoordinateSystem::toXyz(inTheta, inPhi);

    // Get the outgoing direction.
    float outTheta = toRadian(outgoingPolarAngle);
    float outPhi   = toRadian(outgoingAzimuthalAngle);
    Vec3 outDir = SphericalCoordinateSystem::toXyz(outTheta, outPhi);

    if (attributeTypeName == ValueName) {
        // Sample and display the value of BRDF/BTDF.
        Spectrum value = bsdf->getSpectrum(inDir, outDir);
        if (value.size() > 0) {
            std::cout << "Value: " << value.format(LB_EIGEN_IO_FMT) << std::endl;
        }
    }
    else if (attributeTypeName == ReflectanceName) {
        Vec3 sumXyz(0.0, 0.0, 0.0);
        int numComponents = 0;

        if (bsdf) {
            // Compute and display reflectance.
            if (std::shared_ptr<const Brdf> brdf = bsdf->getBrdf()) {
                std::cout << "[BRDF]" << std::endl;
                Vec3 xyz;
                showReflectance(brdf.get(), inDir, &xyz);
                sumXyz += xyz;
                ++numComponents;
            }

            // Compute and display transmittance.
            if (std::shared_ptr<const Btdf> btdf = bsdf->getBtdf()) {
                std::cout << "[BTDF]" << std::endl;
                Vec3 xyz;
                showReflectance(btdf->getBrdf().get(), inDir, &xyz);
                sumXyz += xyz;
                ++numComponents;
            }
        }

        // Compute and display reflectance.
        if (sr) {
            std::cout << "[Specular reflectance]" << std::endl;
            Vec3 xyz;
            showReflectance(sr.get(), inDir, &xyz);
            sumXyz += xyz;
            ++numComponents;
        }

        // Compute and display transmittance.
        if (st) {
            std::cout << "[Specular transmittance]" << std::endl;
            Vec3 xyz;
            showReflectance(st.get(), inDir, &xyz);
            sumXyz += xyz;
            ++numComponents;
        }

        if (numComponents >= 2) {
            std::cout << "[Total]" << std::endl;
            Vec3 sumRgb = xyzToSrgb(sumXyz).cwiseMax(0.0);
            showReflectance(sumXyz, sumRgb);
        }
    }

    return 0;
}
