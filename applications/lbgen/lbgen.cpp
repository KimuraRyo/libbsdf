// =================================================================== //
// Copyright (C) 2018-2019 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#include <iostream>

#include <libbsdf/Brdf/Processor.h>
#include <libbsdf/Brdf/SpecularCoordinatesBrdf.h>

#include <libbsdf/Common/Log.h>
#include <libbsdf/Common/Version.h>

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
    Log::setNotificationLevel(Log::Level::ERROR_MSG);

    ArgumentParser ap(argc, argv);

    using std::cout;
    using std::cerr;
    using std::endl;

    if (ap.read("-h") || ap.read("--help") || ap.getTokens().empty()) {
        cout << "Usage: lbgen [options ...] bsdf_model out_file" << endl;
        cout << endl;
        cout << "lbgen generates BRDF/BTDF data and save an Integra BSDF file." << endl;
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
        cout << "  -numIncomingPolarAngles      set the division number of incoming polar angles (default: 90)" << endl;
        cout << "  -numSpecularPolarAngles      set the division number of specular polar angles (default: 90)" << endl;
        cout << "  -numSpecularAzimuthalAngles  set the division number of specular azimuthal angles (default: 72)" << endl;
        cout << "  -conservationOfEnergy        fix BSDF/BRDF/BTDF if the sum of reflectances and transmittances exceed one" << endl;
#ifdef _OPENMP
        cout << "  -numThreads                  set the number of threads used by parallel processing" << endl;
#endif
        return 0;
    }

    const std::string version("1.0.8");

    if (ap.read("-v") || ap.read("--version")) {
        cout << "Version: lbgen " << version << " (libbsdf-" << getVersion() << ")" << endl;
        return 0;
    }

    const std::string GgxName                       = "ggx";
    const std::string MultipleScatteringSmithName   = "multiple-scattering-smith";
    const std::string LambertianName                = "lambertian";

    if (ap.read("-l") || ap.read("--list")) {
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
        return 0;
    }

    if (ap.read("-r") || ap.read("--reference")) {
        cout << "  " << GgxName << ": " << Ggx::getReference() << endl;
        cout << "  " << MultipleScatteringSmithName << ": " << MultipleScatteringSmith::getReference() << endl;
        return 0;
    }

    int numIncomingPolarAngles = 90;
    if (ap.read("-numIncomingPolarAngles", &numIncomingPolarAngles) == ArgumentParser::ERROR) {
        return 1;
    }
    else {
        numIncomingPolarAngles = utility::clampParameter("numIncomingPolarAngles",
                                                         numIncomingPolarAngles + 1, 2, 3600);
    }

    int numSpecularPolarAngles = 90;
    if (ap.read("-numSpecularPolarAngles", &numSpecularPolarAngles) == ArgumentParser::ERROR) {
        return 1;
    }
    else {
        numSpecularPolarAngles = utility::clampParameter("numSpecularPolarAngles",
                                                         numSpecularPolarAngles + 1, 2, 3600);
    }

    int numSpecularAzimuthalAngles = 72;
    if (ap.read("-numSpecularAzimuthalAngles", &numSpecularAzimuthalAngles) == ArgumentParser::ERROR) {
        return 1;
    }
    else {
        numSpecularAzimuthalAngles = utility::clampParameter("numSpecularAzimuthalAngles",
                                                             numSpecularAzimuthalAngles + 1, 2, 3600);
    }

    bool conservationOfEnergyUsed = false;
    if (ap.read("-conservationOfEnergy")) {
        conservationOfEnergyUsed = true;
    }

#ifdef _OPENMP
    int numThreads;
    if (ap.read("-numThreads", &numThreads) == ArgumentParser::ERROR) {
        return 1;
    }
    else {
        omp_set_num_threads(numThreads);
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
        cerr << "Invalid value (n): " << n << endl;
        return 1;
    }

    float k = 0.0f;
    if (ap.read("-k", &k) == ArgumentParser::ERROR) {
        return 1;
    }
    else if (k < 0.0f) {
        cerr << "Invalid value (k): " << k << endl;
        return 1;
    }

    int numIterations = 10;
    if (ap.read("-numIterations", &numIterations) == ArgumentParser::ERROR) {
        return 1;
    }
    else {
        numIterations = utility::clampParameter("numIterations", numIterations, 1, 10000);
    }

    if (!ap.validate(2)) return 1;

    std::string modelName = ap.getTokens().at(0);
    std::string fileName  = ap.getTokens().at(1);

    ReflectanceModel* model;

    // Create a BRDF/BTDF model.
    if (modelName == GgxName) {
        roughness = utility::clampParameter("roughness", roughness, 0.01f, 1.0f);
        model = new Ggx(Vec3(1.0, 1.0, 1.0), roughness, n, k);
    }
    else if (modelName == MultipleScatteringSmithName) {
        roughness = utility::clampParameter("roughness", roughness, 0.01f, 1.0f);

        MultipleScatteringSmith::MaterialType matType = (k == 0.0f) ? MultipleScatteringSmith::DIELECTRIC_MATERIAL
                                                                    : MultipleScatteringSmith::CONDUCTOR_MATERIAL;

        model = new MultipleScatteringSmith(Vec3(1.0, 1.0, 1.0), roughness, roughness, n,
                                            static_cast<int>(matType),
                                            static_cast<int>(MultipleScatteringSmith::GAUSSIAN_HEIGHT),
                                            static_cast<int>(MultipleScatteringSmith::BECKMANN_SLOPE),
                                            numIterations);
    }
    else if (modelName == LambertianName) {
        n = 1.0f;
        model = new Lambertian(Vec3(1.0, 1.0, 1.0));
    }
    else {
        cerr << "Invalid model name: " << modelName << endl;
        return 1;
    }

    std::string comments("Software: lbgen-" + version);
    comments += "\n;; Arguments:";
    for (int i = 1; i < argc; ++i) {
        comments += " " + std::string(argv[i]);
    }

    // Create BRDF/BTDF and write file.
    if (reader_utility::hasSuffix(fileName, ".ddr")) {
        SpecularCoordinatesBrdf* brdf = createBrdf(*model,
                                                   1.0f,
                                                   numIncomingPolarAngles,
                                                   numSpecularPolarAngles,
                                                   numSpecularAzimuthalAngles,
                                                   BRDF_DATA);

        if (conservationOfEnergyUsed) {
            fixEnergyConservation(brdf);
        }

        DdrWriter::write(fileName, *brdf, comments);

        delete brdf;

        cout << "Generated: " << fileName << endl;
    }
    else if (reader_utility::hasSuffix(fileName, ".ddt")) {
        SpecularCoordinatesBrdf* btdf = createBrdf(*model,
                                                   n,
                                                   numIncomingPolarAngles,
                                                   numSpecularPolarAngles,
                                                   numSpecularAzimuthalAngles,
                                                   BTDF_DATA);

        if (conservationOfEnergyUsed) {
            fixEnergyConservation(btdf);
        }

        DdrWriter::write(fileName, *btdf, comments);

        delete btdf;

        cout << "Generated: " << fileName << endl;
    }
    else {
        SpecularCoordinatesBrdf* brdf = createBrdf(*model,
                                                   1.0f,
                                                   numIncomingPolarAngles,
                                                   numSpecularPolarAngles,
                                                   numSpecularAzimuthalAngles,
                                                   BRDF_DATA);

        SpecularCoordinatesBrdf* btdf = createBrdf(*model,
                                                   n,
                                                   numIncomingPolarAngles,
                                                   numSpecularPolarAngles,
                                                   numSpecularAzimuthalAngles,
                                                   BTDF_DATA);

        if (conservationOfEnergyUsed) {
            fixEnergyConservation(brdf, btdf);
        }

        DdrWriter::write(fileName + ".ddr", *brdf, comments);
        DdrWriter::write(fileName + ".ddt", *btdf, comments);

        cout << "Generated: " << fileName + ".ddr" << endl;
        cout << "Generated: " << fileName + ".ddt" << endl;

        delete brdf;
        delete btdf;
    }

    delete model;

    return 0;
}
