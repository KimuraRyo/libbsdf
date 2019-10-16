// =================================================================== //
// Copyright (C) 2018-2019 Kimura Ryo                                  //
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

#include <ArgumentParser.h>
#include <Utility.h>

using namespace lb;

/*
 * BRDF/BTDF generator
 */

const std::string APP_NAME("lbgen");
const std::string APP_VERSION("1.0.10");

const std::string GgxName                       = "ggx";
const std::string MultipleScatteringSmithName   = "multiple-scattering-smith";
const std::string LambertianName                = "lambertian";

void showHelp()
{
    using std::cout;
    using std::endl;

    cout << "Usage: lbgen [options ...] bsdf_model out_file" << endl;
    cout << endl;
    cout << "lbgen generates BRDF/BTDF data and saves an Integra BSDF file." << endl;
    cout << endl;
    cout << "Positional Arguments:" << endl;
    cout << "  bsdf_model  Analytic BSDF model" << endl;
    cout << "  out_file    Name of generated file";
    cout << " (\".ddr\" or \".ddt\" is acceptable as a suffix for BRDF or BTDF.";
    cout << " Otherwise, \".ddr\" and \".ddt\" files are generated.)" << endl;
    cout << endl;
    cout << "Options:" << endl;
    cout << "  -h, --help                   show this help message and exit" << endl;
    cout << "  -v, --version                show program's version number and exit" << endl;
    cout << "  -l, --list                   show acceptable models and parameters of BSDF and exit" << endl;
    cout << "  -r, --reference              show the references of BSDF models and exit" << endl;
    cout << "  -numIncomingPolarAngles      set the division number of incoming polar angles (default: 18)" << endl;
    cout << "  -numSpecularPolarAngles      set the division number of specular polar angles (default: 360)" << endl;
    cout << "  -numSpecularAzimuthalAngles  set the division number of specular azimuthal angles (default: 72)" << endl;
    cout << "  -conservationOfEnergy        fix BSDF/BRDF/BTDF if the sum of reflectances and transmittances exceed one" << endl;
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

    int numIncomingPolarAngles = 18;
    if (ap.read("-numIncomingPolarAngles", &numIncomingPolarAngles) == ArgumentParser::ERROR) {
        return 1;
    }
    else {
        numIncomingPolarAngles = app_utility::clampParameter("numIncomingPolarAngles",
                                                             numIncomingPolarAngles + 1, 2, 3600);
    }

    int numSpecularPolarAngles = 360;
    if (ap.read("-numSpecularPolarAngles", &numSpecularPolarAngles) == ArgumentParser::ERROR) {
        return 1;
    }
    else {
        numSpecularPolarAngles = app_utility::clampParameter("numSpecularPolarAngles",
                                                             numSpecularPolarAngles + 1, 2, 3600);
    }

    int numSpecularAzimuthalAngles = 72;
    if (ap.read("-numSpecularAzimuthalAngles", &numSpecularAzimuthalAngles) == ArgumentParser::ERROR) {
        return 1;
    }
    else {
        numSpecularAzimuthalAngles = app_utility::clampParameter("numSpecularAzimuthalAngles",
                                                                 numSpecularAzimuthalAngles + 1, 2, 3600);
    }

    bool conservationOfEnergyUsed = false;
    if (ap.read("-conservationOfEnergy")) {
        conservationOfEnergyUsed = true;
    }

#ifdef _OPENMP
    int numThreads;
    ArgumentParser::ResultType result_numThreads = ap.read("-numThreads", &numThreads);
    if (result_numThreads == ArgumentParser::OK) {
        omp_set_num_threads(numThreads);
    }
    else if (result_numThreads == ArgumentParser::ERROR) {
        return 1;
    }
#endif

    float roughness = 0.3f;
    if (ap.read("-roughness", &roughness) == ArgumentParser::ERROR) {
        return 1;
    }

    float n = 1.5f;
    if (ap.read("-n", &n) == ArgumentParser::ERROR) {
        return 1;
    }
    else if (n <= 0.0f) {
        std::cerr << "Invalid value (n): " << n << std::endl;
        return 1;
    }

    float k = 0.0f;
    if (ap.read("-k", &k) == ArgumentParser::ERROR) {
        return 1;
    }
    else if (k < 0.0f) {
        std::cerr << "Invalid value (k): " << k << std::endl;
        return 1;
    }

    int numIterations = 10;
    if (ap.read("-numIterations", &numIterations) == ArgumentParser::ERROR) {
        return 1;
    }
    else {
        numIterations = app_utility::clampParameter("numIterations", numIterations, 1, 10000);
    }

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

    // Create BRDFs/BTDFs and save files.
    std::string comments = app_utility::createComments(argc, argv, APP_NAME, APP_VERSION);
    if (reader_utility::hasSuffix(fileName, ".ddr")) {
        std::unique_ptr<SpecularCoordinatesBrdf> brdf(createBrdf(*model,
                                                                 1.0f,
                                                                 numIncomingPolarAngles,
                                                                 numSpecularPolarAngles,
                                                                 numSpecularAzimuthalAngles,
                                                                 BRDF_DATA));

        if (conservationOfEnergyUsed) {
            fixEnergyConservation(brdf.get());
        }

        if (DdrWriter::write(fileName, *brdf, comments)) {
            std::cout << "Saved: " << fileName << std::endl;
        }
    }
    else if (reader_utility::hasSuffix(fileName, ".ddt")) {
        std::unique_ptr<SpecularCoordinatesBrdf> btdf(createBrdf(*model,
                                                                 n,
                                                                 numIncomingPolarAngles,
                                                                 numSpecularPolarAngles,
                                                                 numSpecularAzimuthalAngles,
                                                                 BTDF_DATA));

        if (conservationOfEnergyUsed) {
            fixEnergyConservation(btdf.get());
        }

        if (DdrWriter::write(fileName, *btdf, comments)) {
            std::cout << "Saved: " << fileName << std::endl;
        }
    }
    else {
        std::unique_ptr<SpecularCoordinatesBrdf> brdf(createBrdf(*model,
                                                                 1.0f,
                                                                 numIncomingPolarAngles,
                                                                 numSpecularPolarAngles,
                                                                 numSpecularAzimuthalAngles,
                                                                 BRDF_DATA));

        std::unique_ptr<SpecularCoordinatesBrdf> btdf(createBrdf(*model,
                                                                 n,
                                                                 numIncomingPolarAngles,
                                                                 numSpecularPolarAngles,
                                                                 numSpecularAzimuthalAngles,
                                                                 BTDF_DATA));

        if (conservationOfEnergyUsed) {
            fixEnergyConservation(brdf.get(), btdf.get());
        }

        if (DdrWriter::write(fileName + ".ddr", *brdf, comments)) {
            std::cout << "Saved: " << fileName + ".ddr" << std::endl;
        }

        if (DdrWriter::write(fileName + ".ddt", *btdf, comments)) {
            std::cout << "Saved: " << fileName + ".ddt" << std::endl;
        }
    }

    return 0;
}
