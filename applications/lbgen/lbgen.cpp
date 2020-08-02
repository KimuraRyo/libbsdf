// =================================================================== //
// Copyright (C) 2018-2020 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#include <iostream>
#include <memory>

#include <libbsdf/Brdf/Processor.h>

#include <libbsdf/Reader/ReaderUtility.h>

#include <libbsdf/ReflectanceModel/GGX.h>
#include <libbsdf/ReflectanceModel/Lambertian.h>
#include <libbsdf/ReflectanceModel/MultipleScatteringSmith.h>
#include <libbsdf/ReflectanceModel/ReflectanceModelUtility.h>

#include <libbsdf/Writer/DdrWriter.h>
#include <libbsdf/Writer/SsddWriter.h>

#include <ArgumentParser.h>
#include <Utility.h>

using namespace lb;

/*
 * BSDF generator
 */

const std::string APP_NAME("lbgen");
const std::string APP_VERSION("1.0.11");

const std::string GgxName                       = "ggx";
const std::string MultipleScatteringSmithName   = "multiple-scattering-smith";
const std::string LambertianName                = "lambertian";

// Parameters
int numIncomingPolarAngles = 18;
int numSpecularPolarAngles = 360;
int numSpecularAzimuthalAngles = 72;
bool conservationOfEnergyUsed = false;
float roughness = 0.3f;
float n = 1.5f;
float k = 0.0f;
int numIterations = 10;

void showHelp()
{
    using std::cout;
    using std::endl;

    cout << "Usage: lbgen [options ...] bsdf_model out_file" << endl;
    cout << endl;
    cout << "lbgen generates BSDF data and saves the file(s)." << endl;
    cout << endl;
    cout << "Positional Arguments:" << endl;
    cout << "  bsdf_model   Analytic BSDF model" << endl;
    cout << "  out_file     Name of generated file";
    cout << "                   Valid formats:" << endl;
    cout << "                       Surface Scattering Distribution Data (\".ssdd\")" << endl;
    cout << "                       Integra Diffuse Distribution (\".ddr, .ddt\")" << endl;
    cout << "                           If there is no file extension, \".ddr\" and \".ddt\" files are generated." << endl;
    cout << endl;
    cout << "Options:" << endl;
    cout << "  -h, --help                   show this help message and exit" << endl;
    cout << "  -v, --version                show program's version number and exit" << endl;
    cout << "  -l, --list                   show acceptable models and parameters of BSDF and exit" << endl;
    cout << "  -r, --reference              show the references of BSDF models and exit" << endl;
    cout << "  -numIncomingPolarAngles      set the division number of incoming polar angles (default: 18)" << endl;
    cout << "  -numSpecularPolarAngles      set the division number of specular polar angles (default: 360)" << endl;
    cout << "  -numSpecularAzimuthalAngles  set the division number of specular azimuthal angles (default: 72)" << endl;
    cout << "  -conservationOfEnergy        fix BSDF/BRDF/BTDF if the sum of reflectance and transmittance exceed one" << endl;
#ifdef _OPENMP
    cout << "  -numThreads                  set the number of threads used by parallel processing" << endl;
#endif
}

void showList()
{
    using std::cout;
    using std::endl;

    cout << "Acceptable models:" << endl;
    cout << "  " << GgxName << endl;
    cout << "      Valid options:" << endl;
    cout << "          -roughness  set roughness of surface (default: 0.3, range: [0.01, 1.0])" << endl;
    cout << "          -n          set refractive index (default: 1.5)" << endl;
    cout << "          -k          set extinction coefficient (default: 0.0)" << endl;
    cout << "  " << MultipleScatteringSmithName << endl;
    cout << "      Valid options:" << endl;
    cout << "          -roughness      set roughness of surface (default: 0.3, range: [0.01, 1.0])" << endl;
    cout << "          -n              set refractive index (default: 1.5)" << endl;
    cout << "          -numIterations  set the number of sampling iterations (default: 10)" << endl;
    cout << "  " << LambertianName << endl;
    cout << "      Valid options: none" << endl;
}

void showReference()
{
    using std::cout;
    using std::endl;

    cout << "  " << GgxName << ": " << Ggx::getReference() << endl;
    cout << "  " << MultipleScatteringSmithName << ": " << MultipleScatteringSmith::getReference() << endl;
}

bool readOptions(ArgumentParser* ap)
{
    if (ap->read("-numIncomingPolarAngles", &numIncomingPolarAngles) == ArgumentParser::ERROR) {
        return false;
    }
    else {
        numIncomingPolarAngles = app_utility::clampParameter("numIncomingPolarAngles",
                                                             numIncomingPolarAngles + 1, 2, 3600);
    }

    if (ap->read("-numSpecularPolarAngles", &numSpecularPolarAngles) == ArgumentParser::ERROR) {
        return false;
    }
    else {
        numSpecularPolarAngles = app_utility::clampParameter("numSpecularPolarAngles",
                                                             numSpecularPolarAngles + 1, 2, 3600);
    }

    if (ap->read("-numSpecularAzimuthalAngles", &numSpecularAzimuthalAngles) == ArgumentParser::ERROR) {
        return false;
    }
    else {
        numSpecularAzimuthalAngles = app_utility::clampParameter("numSpecularAzimuthalAngles",
                                                                 numSpecularAzimuthalAngles + 1, 2, 3600);
    }

    if (ap->read("-conservationOfEnergy")) {
        conservationOfEnergyUsed = true;
    }

#ifdef _OPENMP
    int numThreads;
    ArgumentParser::ResultType result_numThreads = ap->read("-numThreads", &numThreads);
    if (result_numThreads == ArgumentParser::OK) {
        omp_set_num_threads(numThreads);
    }
    else if (result_numThreads == ArgumentParser::ERROR) {
        return false;
    }
#endif

    if (ap->read("-roughness", &roughness) == ArgumentParser::ERROR) {
        return false;
    }

    if (ap->read("-n", &n) == ArgumentParser::ERROR) {
        return false;
    }
    else if (n <= 0.0f) {
        std::cerr << "Invalid value (n): " << n << std::endl;
        return false;
    }

    if (ap->read("-k", &k) == ArgumentParser::ERROR) {
        return false;
    }
    else if (k < 0.0f) {
        std::cerr << "Invalid value (k): " << k << std::endl;
        return false;
    }

    if (ap->read("-numIterations", &numIterations) == ArgumentParser::ERROR) {
        return false;
    }
    else {
        numIterations = app_utility::clampParameter("numIterations", numIterations, 1, 10000);
    }

    return true;
}

SpecularCoordinatesBrdf* createBrdf(const ReflectanceModel& model,
                                    float                   n,
                                    int                     numIncomingPolarAngles,
                                    int                     numSpecularPolarAngles,
                                    int                     numSpecularAzimuthalAngles,
                                    DataType                dataType)
{
    n = std::max(n, 1.0f);

    SpecularCoordinatesBrdf* brdf = new SpecularCoordinatesBrdf(numIncomingPolarAngles,
                                                                1,
                                                                numSpecularPolarAngles,
                                                                numSpecularAzimuthalAngles,
                                                                2.0f,
                                                                MONOCHROMATIC_MODEL, 1, n);

    reflectance_model_utility::setupTabularBrdf(model, brdf, dataType);
    brdf->setSourceType(GENERATED_SOURCE);

    return brdf;
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

    if (ap.read("-r") || ap.read("--reference")) {
        showReference();
        return 0;
    }

    if (!readOptions(&ap)) return 1;

    if (!ap.validateNumTokens(2)) return 1;

    std::string modelName = ap.getTokens().at(0);
    std::string fileName  = ap.getTokens().at(1);

    // Create a BRDF/BTDF model.
    std::unique_ptr<ReflectanceModel> model;
    const Vec3 white(1.0, 1.0, 1.0);
    if (modelName == GgxName) {
        roughness = app_utility::clampParameter("roughness", roughness, 0.01f, 1.0f);
        model.reset(new Ggx(white, roughness, n, k));
    }
    else if (modelName == MultipleScatteringSmithName) {
        roughness = app_utility::clampParameter("roughness", roughness, 0.01f, 1.0f);

        MultipleScatteringSmith::MaterialType matType = (k == 0.0f) ? MultipleScatteringSmith::DIELECTRIC_MATERIAL
                                                                    : MultipleScatteringSmith::CONDUCTOR_MATERIAL;

        model.reset(new MultipleScatteringSmith(white, roughness, roughness, n,
                                                static_cast<int>(matType),
                                                static_cast<int>(MultipleScatteringSmith::GAUSSIAN_HEIGHT),
                                                static_cast<int>(MultipleScatteringSmith::BECKMANN_SLOPE),
                                                numIterations));
    }
    else if (modelName == LambertianName) {
        n = 1.0f;
        model.reset(new Lambertian(white));
    }
    else {
        std::cerr << "Invalid model name: " << modelName << std::endl;
        return 1;
    }

    FileType fileType;
    if (reader_utility::hasSuffix(fileName, ".ssdd")) {
        fileType = SSDD_FILE;
    }
    else if (reader_utility::hasSuffix(fileName, ".ddr")) {
        fileType = INTEGRA_DDR_FILE;
    }
    else if (reader_utility::hasSuffix(fileName, ".ddt")) {
        fileType = INTEGRA_DDT_FILE;
    }

    std::unique_ptr<SpecularCoordinatesBrdf> brdf, btdfData;

    // Create a BRDF.
    if (fileType != INTEGRA_DDT_FILE) {
        brdf.reset(createBrdf(*model,
                              1.0f,
                              numIncomingPolarAngles,
                              numSpecularPolarAngles,
                              numSpecularAzimuthalAngles,
                              BRDF_DATA));
    }

    // Create BTDF data in lb::SpecularCoordinatesBrdf.
    if (fileType != INTEGRA_DDR_FILE) {
        btdfData.reset(createBrdf(*model,
                                  n,
                                  numIncomingPolarAngles,
                                  numSpecularPolarAngles,
                                  numSpecularAzimuthalAngles,
                                  BTDF_DATA));
    }

    // Save files.
    std::string comments = app_utility::createComments(argc, argv, APP_NAME, APP_VERSION);
    if (fileType == INTEGRA_DDR_FILE) {
        if (conservationOfEnergyUsed) {
            fixEnergyConservation(brdf.get());
        }

        if (DdrWriter::write(fileName, *brdf, comments)) {
            std::cout << "Saved: " << fileName << std::endl;
        }
    }
    else if (fileType == INTEGRA_DDT_FILE) {
        if (conservationOfEnergyUsed) {
            fixEnergyConservation(btdfData.get());
        }

        if (DdrWriter::write(fileName, *btdfData, comments)) {
            std::cout << "Saved: " << fileName << std::endl;
        }
    }
    else {
        if (conservationOfEnergyUsed) {
            fixEnergyConservation(brdf.get(), btdfData.get());
        }

        if (fileType == SSDD_FILE) {
            std::unique_ptr<Btdf> btdf(new Btdf(std::move(btdfData)));
            std::unique_ptr<Bsdf> bsdf(new Bsdf(std::move(brdf), std::move(btdf)));
            std::unique_ptr<Material> material(new Material(std::move(bsdf)));

            if (SsddWriter::write(fileName, *material, SsddWriter::DataFormat::ASCII_DATA, comments)) {
                std::cout << "Saved: " << fileName << std::endl;
            }
        }
        else {
            if (DdrWriter::write(fileName + ".ddr", *brdf, comments)) {
                std::cout << "Saved: " << fileName + ".ddr" << std::endl;
            }

            if (DdrWriter::write(fileName + ".ddt", *btdfData, comments)) {
                std::cout << "Saved: " << fileName + ".ddt" << std::endl;
            }
        }
    }

    return 0;
}
