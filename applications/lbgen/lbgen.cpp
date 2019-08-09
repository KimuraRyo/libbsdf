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

#include <libbsdf/Common/Version.h>

#include <libbsdf/Reader/ReaderUtility.h>

#include <libbsdf/ReflectanceModel/GGX.h>
#include <libbsdf/ReflectanceModel/Lambertian.h>
#include <libbsdf/ReflectanceModel/MultipleScatteringSmith.h>
#include <libbsdf/ReflectanceModel/ReflectanceModelUtility.h>

#include <libbsdf/Writer/DdrWriter.h>

#include <ArgumentParser.h>

using namespace lb;

/*
 * BRDF/BTDF generator
 */

template <typename T>
T clampParameter(const std::string& name, T val, T minVal, T maxVal)
{
    if (val < minVal || val > maxVal) {
        val = clamp(val, minVal, maxVal);
        std::cout << "Warning: \"" + name + "\" is clamped to " << val << "." << std::endl;
    }

    return val;
}

SpecularCoordinatesBrdf* createBrdf(const ReflectanceModel& model,
                                    float                   n,
                                    int                     numIncomingPolarAngles,
                                    int                     numSpecularPolarAngles,
                                    int                     numSpecularAzimuthalAngles,
                                    DataType                dataType)
{
    std::cout.setstate(std::ios_base::failbit);
    SpecularCoordinatesBrdf* brdf = new SpecularCoordinatesBrdf(numIncomingPolarAngles,
                                                                1,
                                                                numSpecularPolarAngles,
                                                                numSpecularAzimuthalAngles,
                                                                2.0f,
                                                                MONOCHROMATIC_MODEL, 1, n);
    std::cout.clear();

    reflectance_model_utility::setupTabularBrdf(model, brdf, dataType);
    brdf->setSourceType(GENERATED_SOURCE);

    return brdf;
}

int main(int argc, char** argv)
{
    ArgumentParser ap(argc, argv);

    if (ap.read("-h") || ap.read("--help") || ap.getTokens().empty()) {
        std::cout << "Usage: lbgen [options ...] bsdf_model out_file" << std::endl;
        std::cout << std::endl;
        std::cout << "lbgen generates BRDF/BTDF data and save an Integra BSDF file." << std::endl;
        std::cout << std::endl;
        std::cout << "Positional Arguments:" << std::endl;
        std::cout << "  bsdf_model  Analytic BSDF model" << std::endl;
        std::cout << "  out_file    Name of generated file";
        std::cout << " (\".ddr\" or \".ddt\" is acceptable as a suffix for BRDF or BTDF.";
        std::cout << " Otherwise, \".ddr\" and \".ddt\" files are generated.)" << std::endl;
        std::cout << std::endl;
        std::cout << "Options:" << std::endl;
        std::cout << "  -h, --help                  show this help message and exit" << std::endl;
        std::cout << "  -v, --version               show program's version number and exit" << std::endl;
        std::cout << "  -l, --list                  show acceptable models and parameters of BSDF and exit" << std::endl;
        std::cout << "  -r, --reference             show the references of BSDF models and exit" << std::endl;
        std::cout << "  -numIncomingPolarAngles     set the division number of incoming polar angles (default: 90)" << std::endl;
        std::cout << "  -numSpecularPolarAngles     set the division number of specular polar angles (default: 90)" << std::endl;
        std::cout << "  -numSpecularAzimuthalAngles set the division number of specular azimuthal angles (default: 72)" << std::endl;
        std::cout << "  -conservationOfEnergy       fix BSDF/BRDF/BTDF if the sum of reflectances and transmittances exceed one" << std::endl;
#ifdef _OPENMP
        std::cout << "  -numThreads                 set the number of threads used by parallel processing" << std::endl;
#endif
        return 0;
    }

    const std::string version("1.0.8");

    if (ap.read("-v") || ap.read("--version")) {
        std::cout << "Version: lbgen " << version << " (libbsdf-" << getVersion() << ")" << std::endl;
        return 0;
    }

    const std::string GgxName                       = "ggx";
    const std::string MultipleScatteringSmithName   = "multiple-scattering-smith";
    const std::string LambertianName                = "lambertian";

    if (ap.read("-l") || ap.read("--list")) {
        std::cout << "Acceptable models:" << std::endl;
        std::cout << "  " << GgxName << std::endl;
        std::cout << "      Valid options:" << std::endl;
        std::cout << "          -roughness  set roughness of surface (default: 0.3, range: [0.01, 1.0])" << std::endl;
        std::cout << "          -n          set refractive index (default: 1.5)" << std::endl;
        std::cout << "          -k          set extinction coefficient (default: 0.0)" << std::endl;
        std::cout << "  " << MultipleScatteringSmithName << std::endl;
        std::cout << "      Valid options:" << std::endl;
        std::cout << "          -roughness      set roughness of surface (default: 0.3, range: [0.01, 1.0])" << std::endl;
        std::cout << "          -n              set refractive index (default: 1.5)" << std::endl;
        std::cout << "          -numIterations  set the number of sampling iterations (default: 10)" << std::endl;
        std::cout << "  " << LambertianName << std::endl;
        std::cout << "      Valid options: none" << std::endl;
        return 0;
    }

    if (ap.read("-r") || ap.read("--reference")) {
        std::cout
            << "  " << GgxName << ": "
            << Ggx::getReference() << std::endl;
        std::cout
            << "  " << MultipleScatteringSmithName << ": "
            << MultipleScatteringSmith::getReference() << std::endl;
        return 0;
    }

    int numIncomingPolarAngles = 90;
    if (ap.read("-numIncomingPolarAngles", &numIncomingPolarAngles) == ArgumentParser::ERROR) {
        return 1;
    }
    else {
        numIncomingPolarAngles = clampParameter("numIncomingPolarAngles",
                                                numIncomingPolarAngles + 1, 2, 3600);
    }

    int numSpecularPolarAngles = 90;
    if (ap.read("-numSpecularPolarAngles", &numSpecularPolarAngles) == ArgumentParser::ERROR) {
        return 1;
    }
    else {
        numSpecularPolarAngles = clampParameter("numSpecularPolarAngles",
                                                numSpecularPolarAngles + 1, 2, 3600);
    }

    int numSpecularAzimuthalAngles = 72;
    if (ap.read("-numSpecularAzimuthalAngles", &numSpecularAzimuthalAngles) == ArgumentParser::ERROR) {
        return 1;
    }
    else {
        numSpecularAzimuthalAngles = clampParameter("numSpecularAzimuthalAngles",
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
        numIterations = clampParameter("numIterations", numIterations, 1, 10000);
    }

    if (!ap.validate(2)) return 1;

    std::string modelName = ap.getTokens().at(0);
    std::string fileName  = ap.getTokens().at(1);

    ReflectanceModel* model;

    // Create a BRDF/BTDF model.
    if (modelName == GgxName) {
        roughness = clampParameter("roughness", roughness, 0.01f, 1.0f);
        model = new Ggx(Vec3(1.0, 1.0, 1.0), roughness, n, k);
    }
    else if (modelName == MultipleScatteringSmithName) {
        roughness = clampParameter("roughness", roughness, 0.01f, 1.0f);

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
        std::cerr << "Invalid model name: " << modelName << std::endl;
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

        std::cout << "Generated: " << fileName << std::endl;
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

        std::cout << "Generated: " << fileName << std::endl;
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
            std::cout.setstate(std::ios_base::failbit);
            fixEnergyConservation(brdf, btdf);
            std::cout.clear();
        }

        DdrWriter::write(fileName + ".ddr", *brdf, comments);
        DdrWriter::write(fileName + ".ddt", *btdf, comments);

        std::cout << "Generated: " << fileName + ".ddr" << std::endl;
        std::cout << "Generated: " << fileName + ".ddt" << std::endl;

        delete brdf;
        delete btdf;
    }

    delete model;

    return 0;
}
