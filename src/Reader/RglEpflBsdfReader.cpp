// =================================================================== //
// Copyright (C) 2026 Kimura Ryo                                       //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#include <libbsdf/Reader/RglEpflBsdfReader.h>

#define POWITACQ_IMPLEMENTATION
#include <src/ThirdParty/rgl-epfl_brdf-loader/powitacq.h>
#include <src/ThirdParty/rgl-epfl_brdf-loader/powitacq_rgb.h>

#include <libbsdf/Brdf/Processor.h>
#include <libbsdf/Common/Log.h>

using namespace lb;

SpecularCoordinatesBrdf* convertPowitacqToBrdf(const powitacq::BRDF& powitacqBrdf)
{
    const int                 numInTheta = 18;
    int                       numInPhi = powitacqBrdf.isIsotropic() ? 0 : 72;
    const int                 numSpecTheta = 30;
    const int                 numSpecPhi = 72;
    const powitacq::Spectrum& wavelengths = powitacqBrdf.wavelengths();

    SpecularCoordinatesBrdf* brdf =
        new SpecularCoordinatesBrdf(numInTheta + 1, numInPhi + 1, numSpecTheta + 1, numSpecPhi + 1,
                                    2.0, SPECTRAL_MODEL, static_cast<int>(wavelengths.size()));

    SampleSet* ss = brdf->getSampleSet();

    for (int i = 0; i < static_cast<int>(wavelengths.size()); ++i) {
        ss->setWavelength(i, wavelengths[i]);
    }

    for (int i0 = 0; i0 < ss->getNumAngles0(); ++i0) {
    for (int i1 = 0; i1 < ss->getNumAngles1(); ++i1) {
    for (int i2 = 0; i2 < ss->getNumAngles2(); ++i2) {
        Vec3 inDir, outDir;
        powitacq::Spectrum powitacqSp;
        #pragma omp parallel for private(inDir, outDir, powitacqSp)
        for (int i3 = 0; i3 < ss->getNumAngles3(); ++i3) {
            brdf->getInOutDirection(i0, i1, i2, i3, &inDir, &outDir);

            powitacqSp = powitacqBrdf.eval(powitacq::Vector3f(inDir.x(), inDir.y(), inDir.z()),
                                           powitacq::Vector3f(outDir.x(), outDir.y(), outDir.z()));

            Eigen::ArrayXf sp = Eigen::Map<Eigen::ArrayXf>(&powitacqSp[0], powitacqSp.size());
            for (auto& spVal : sp) {
                if (std::isnan(spVal))
                    spVal = 0.0f;
            }
            brdf->setSpectrum(i0, i1, i2, i3, sp);
        }
    }}}

    return brdf;
}

SpecularCoordinatesBrdf* convertPowitacqRgbToBrdf(const powitacq_rgb::BRDF& powitacqRgbBrdf)
{
    const int numInTheta = 18;
    int       numInPhi = powitacqRgbBrdf.isIsotropic() ? 0 : 72;
    const int numSpecTheta = 30;
    const int numSpecPhi = 72;

    SpecularCoordinatesBrdf* brdf = new SpecularCoordinatesBrdf(
        numInTheta + 1, numInPhi + 1, numSpecTheta + 1, numSpecPhi + 1, 2.0, RGB_MODEL);

    SampleSet* ss = brdf->getSampleSet();
    for (int i0 = 0; i0 < ss->getNumAngles0(); ++i0) {
    for (int i1 = 0; i1 < ss->getNumAngles1(); ++i1) {
    for (int i2 = 0; i2 < ss->getNumAngles2(); ++i2) {
        Vec3 inDir, outDir;
        powitacq_rgb::Vector3f powitacqRgb;
        #pragma omp parallel for private(inDir, outDir, powitacqRgb)
        for (int i3 = 0; i3 < ss->getNumAngles3(); ++i3) {
            brdf->getInOutDirection(i0, i1, i2, i3, &inDir, &outDir);

            powitacqRgb =
                powitacqRgbBrdf.eval(powitacq_rgb::Vector3f(inDir.x(), inDir.y(), inDir.z()),
                                     powitacq_rgb::Vector3f(outDir.x(), outDir.y(), outDir.z()));

            Vec3f rgb(powitacqRgb[0], powitacqRgb[1], powitacqRgb[2]);
            brdf->setSpectrum(i0, i1, i2, i3, rgb);
        }
    }}}

    return brdf;
}

SpecularCoordinatesBrdf* RglEpflBsdfReader::read(const std::string& fileName)
{
    SpecularCoordinatesBrdf* brdf = nullptr;

    try {
        powitacq::BRDF powitacqBrdf(fileName);
        brdf = convertPowitacqToBrdf(powitacqBrdf);
    }
    catch (const std::runtime_error& e) {
        // If spectral data can not be read, RGB data is read.
        if (e.what() != std::string("Tensor: Unable to find field spectra")) {
            lbError << "[RglEpflBsdfReader::read] Failed to read: " << fileName << ", " << e.what();
            return nullptr;
        }
    }
    catch (const std::exception& e) {
        lbError << "[RglEpflBsdfReader::read] Failed to read: " << fileName << ", " << e.what();
        return nullptr;
    }

    if (!brdf) {
        try {
            powitacq_rgb::BRDF powitacqRgbBrdf(fileName);
            brdf = convertPowitacqRgbToBrdf(powitacqRgbBrdf);
        }
        catch (const std::exception& e) {
            lbError << "[RglEpflBsdfReader::read] Failed to read: " << fileName << ", " << e.what();
            return nullptr;
        }
    }

    // Convert CCBRDF values to BRDF.
    divideByCosineOutTheta(brdf, 0.001);

    // Extrapolate sample points in areas that have not been measured.
    extrapolateSamplesAlongSpecTheta(brdf, toRadian(82.0), toRadian(84.0));

    fixNegativeSpectra(brdf, false);
    brdf->setSourceType(MEASURED_SOURCE);

    return brdf;
}