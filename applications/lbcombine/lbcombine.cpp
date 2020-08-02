// =================================================================== //
// Copyright (C) 2020 Kimura Ryo                                       //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#include <iostream>
#include <memory>

#include <libbsdf/Reader/DdrReader.h>
#include <libbsdf/Reader/ReaderUtility.h>
#include <libbsdf/Reader/SsddReader.h>
#include <libbsdf/Writer/SsddWriter.h>

#include <ArgumentParser.h>
#include <Utility.h>

using namespace lb;

/*
 * BRDF/BTDF combiner
 */

const std::string APP_NAME("lbcombine");
const std::string APP_VERSION("1.0.0");

// Parameters
std::string brdfFileName;   // BRDF file
std::string btdfFileName;   // BTDF file
std::string srFileName;     // Specular reflectance file
std::string stFileName;     // Specular transmittance file

void showHelp()
{
    using std::cout;
    using std::endl;

    cout << "Usage: lbcombine -data_type1 in_file1 [-data_type2 in_file2] ... out_file" << endl;
    cout << endl;
    cout << "lbcombine combines BRDF, BTDF, specular reflectance, and specular transmittance files into an SSDD file." << endl;
    cout << endl;
    cout << "Positional Arguments:" << endl;
    cout << "  in_file      Name of an input BRDF/BTDF file." << endl;
    cout << "               Valid formats:" << endl;
    cout << "                   Surface Scattering Distribution Data (\".ssdd\")" << endl;
    cout << "                   Integra Diffuse Distribution (\".ddr, .ddt\")" << endl;
    cout << "                   Integra Specular Distribution (\".sdr, .sdt\")" << endl;
    cout << "  out_file     Name of an output BRDF/BTDF file." << endl;
    cout << "               Valid formats:" << endl;
    cout << "                   Surface Scattering Distribution Data (\".ssdd\")" << endl;
    cout << endl;
    cout << "Data Types:" << endl;
    cout << "  -brdf                    append a BRDF to the SSDD file" << endl;
    cout << "  -btdf                    append a BTDF to the SSDD file" << endl;
    cout << "  -specularReflectance     append a specular reflectance to the SSDD file" << endl;
    cout << "  -specularTransmittance   append a specular transmittance to the SSDD file" << endl;
    cout << endl;
    cout << "Options:" << endl;
    cout << "  -h, --help       show this help message and exit" << endl;
    cout << "  -v, --version    show program's version number and exit" << endl;
}

bool readOptions(ArgumentParser* ap)
{
    ap->read("-brdf",                   &brdfFileName);
    ap->read("-btdf",                   &btdfFileName);
    ap->read("-specularReflectance",    &srFileName);
    ap->read("-specularTransmittance",  &stFileName);

    if (brdfFileName.empty() &&
        btdfFileName.empty() &&
        srFileName.empty() &&
        stFileName.empty()) {
        std::cerr << "File names are not found." << std::endl;
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

    if (!readOptions(&ap)) return 1;

    if (!ap.validateNumTokens(1)) return 1;

    std::string outFileName = ap.getTokens().at(0);

    std::shared_ptr<Brdf> brdf;
    std::shared_ptr<Btdf> btdf;
    std::shared_ptr<SampleSet2D> sr;
    std::shared_ptr<SampleSet2D> st;

    // Read a BRDF file.
    if (!brdfFileName.empty()) {
        FileType fileType;
        std::shared_ptr<Material> material = reader_utility::readMaterial(brdfFileName, &fileType);
        if (!material) {
            std::cerr << "Failed to load: " << brdfFileName << std::endl;
            return 1;
        }

        if (std::shared_ptr<Bsdf> bsdf = material->getBsdf()) {
            brdf = bsdf->getBrdf();
        }
    }

    // Read a BTDF file.
    if (!btdfFileName.empty()) {
        FileType fileType;
        std::shared_ptr<Material> material = reader_utility::readMaterial(btdfFileName, &fileType);
        if (!material) {
            std::cerr << "Failed to load: " << btdfFileName << std::endl;
            return 1;
        }

        if (std::shared_ptr<Bsdf> bsdf = material->getBsdf()) {
            btdf = bsdf->getBtdf();
        }
    }

    // Read a specular reflectance file.
    if (!srFileName.empty()) {
        FileType fileType;
        std::shared_ptr<Material> material = reader_utility::readMaterial(srFileName, &fileType);
        if (!material) {
            std::cerr << "Failed to load: " << srFileName << std::endl;
            return 1;
        }

        sr = material->getSpecularReflectances();
    }

    // Read a specular transmittance file.
    if (!stFileName.empty()) {
        FileType fileType;
        std::shared_ptr<Material> material = reader_utility::readMaterial(stFileName, &fileType);
        if (!material) {
            std::cerr << "Failed to load: " << stFileName << std::endl;
            return 1;
        }

        st = material->getSpecularTransmittances();
    }

    std::shared_ptr<Bsdf> outBsdf = std::make_shared<Bsdf>(brdf, btdf);
    std::unique_ptr<Material> outMaterial(new Material(outBsdf, sr, st));

    if (!outMaterial || outMaterial->isEmpty()) {
        std::cerr << "Data is not found." << std::endl;
        return 1;
    }

    // Write an SSDD file.
    std::string comments = app_utility::createComments(argc, argv, APP_NAME, APP_VERSION);
    if (SsddWriter::write(outFileName, *outMaterial, SsddWriter::DataFormat::ASCII_DATA, comments)) {
        std::cout << "Saved: " << outFileName << std::endl;
    }

    return 0;
}
