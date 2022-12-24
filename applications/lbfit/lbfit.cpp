// =================================================================== //
// Copyright (C) 2022 Kimura Ryo                                       //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#include <iostream>
#include <memory>

#include <libbsdf/Brdf/Analyzer.h>
#include <libbsdf/Brdf/Processor.h>
#include <libbsdf/Brdf/SpecularCoordinatesBrdf.h>
#include <libbsdf/Fitter/LambertianFitter.h>
#include <libbsdf/Fitter/SimpleAnisotropicGgxFitter.h>
#include <libbsdf/Fitter/SimpleGgxFitter.h>
#include <libbsdf/Fitter/UnrealEngine4Fitter.h>
#include <libbsdf/Reader/ReaderUtility.h>
#include <libbsdf/ReflectanceModel/ReflectanceModelUtility.h>

#include <ArgumentParser.h>
#include <Utility.h>

using namespace lb;

/*
 * Reflectance model fitter for BRDF
 */

const std::string APP_NAME("lbfit");
const std::string APP_VERSION("1.0.0");

// Parameters
int   numSampling = 100000;
float maxPolarAngle = 75.0f;
bool  fullReport = false;

void showHelp()
{
    using std::cout;
    using std::endl;

    cout << "Usage: lbfit [options ...] file" << endl;
    cout << endl;
    cout << "lbfit fits analytical reflectance models to a loaded BRDF. Fitted parameters and an error are displayed for each reflectance model." << endl;
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
    cout << "  -numSampling     set the number of samples for fitting. If 0, samples in the loaded BRDF are used. (default: 100000)" << endl;
    cout << "  -maxPolarAngle   set the maximum incoming and outgoing polar angle of sample points for fitting. (default: 80.0)" << endl;
    cout << "  -fullReport      show the full report of fitting by Ceres Solver" << endl;
}

bool readOptions(ArgumentParser* ap)
{
    if (ap->read("-numSampling", &numSampling) == ArgumentParser::ERROR) {
        return false;
    }

    if (ap->read("-maxPolarAngle", &maxPolarAngle) == ArgumentParser::ERROR) {
        return false;
    }

    if (ap->read("-fullReport")) {
        fullReport = true;
    }

    return true;
}

void showError(const Brdf& loadedBrdf, const ReflectanceModel& model)
{
    std::unique_ptr<Brdf> loadedRgbBrdf(toSrgb(loadedBrdf));

    std::unique_ptr<Brdf> generatedBrdf(
        new SpecularCoordinatesBrdf(19, 1, 91, 73, 2.0f, RGB_MODEL, 3));
    ReflectanceModelUtility::setupBrdf(model, generatedBrdf.get());

    Spectrum errorSp = computeDifference(*loadedRgbBrdf, *generatedBrdf);
    std::cout << "Error: " << errorSp.sum() / errorSp.size() << std::endl;
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

    std::string fileName = ap.getTokens().at(0);

    FileType fileType;
    DataType dataTyp;
    std::shared_ptr<Brdf> brdf = reader_utility::readBrdf(fileName, &fileType, &dataTyp);
    if (!brdf) {
        std::cerr << "Failed to load: " << fileName << std::endl;
        return 1;
    }

    if (fullReport) {
        Log::setNotificationLevel(Log::Level::INFO_MSG);
    }

    float maxRadian = toRadian(maxPolarAngle);

    UnrealEngine4 ue4 = UnrealEngine4Fitter::estimateParameters(*brdf, numSampling, maxRadian);
    std::cout << "[" << ue4.getName() << "]" << std::endl;
    ReflectanceModelUtility::dumpParametersInfo(ue4);
    showError(*brdf, ue4);

    SimpleGgx simpleGgx = SimpleGgxFitter::estimateParameters(*brdf, numSampling, maxRadian);
    std::cout << "\n[" << simpleGgx.getName() << "]" << std::endl;
    ReflectanceModelUtility::dumpParametersInfo(simpleGgx);
    showError(*brdf, simpleGgx);

    SimpleAnisotropicGgx simpleAnisoGgx =
        SimpleAnisotropicGgxFitter::estimateParameters(*brdf, numSampling, maxRadian);
    std::cout << "\n[" << simpleAnisoGgx.getName() << "]" << std::endl;
    ReflectanceModelUtility::dumpParametersInfo(simpleAnisoGgx);
    showError(*brdf, simpleAnisoGgx);

    Lambertian lambertian = LambertianFitter::estimateParameters(*brdf, numSampling, maxRadian);
    std::cout << "\n[" << lambertian.getName() << "]" << std::endl;
    ReflectanceModelUtility::dumpParametersInfo(lambertian);
    showError(*brdf, lambertian);

    return 0;
}
