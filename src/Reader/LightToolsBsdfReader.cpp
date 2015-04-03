// =================================================================== //
// Copyright (C) 2015 Kimura Ryo                                       //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#include <libbsdf/Brdf/SphericalCoordinatesBrdf.h>
#include <libbsdf/Reader/LightToolsBsdfReader.h>

#include <fstream>
#include <iostream>
#include <set>

using namespace lb;

TwoSidedMaterial* LightToolsBsdfReader::read(const std::string& fileName)
{
    std::ifstream fin(fileName.c_str());
    if (fin.fail()) {
        std::cerr << "[LightToolsBsdfReader::read] Could not open: " << fileName << std::endl;
        return 0;
    }

    std::ios_base::sync_with_stdio(false);

    SymmetryType symmetryType = ASYMMETRICAL;

    ColorModel::Type colorModel = ColorModel::MONOCHROME;

    std::vector<float> outThetaDegrees;
    std::vector<float> outPhiDegrees;

    int numChannels = 1;

    ignoreCommentLines(fin);

    // Read a header.
    std::string headStr;
    while (fin >> headStr) {
        ignoreCommentLines(fin);

        if (headStr.empty()) {
            continue;
        }
        else if (headStr == "Symmetry") {
            std::string typeStr;
            fin >> typeStr;

            if (typeStr == "Asymmetric") {
                symmetryType = ASYMMETRICAL;
            }
        }
        else if (headStr == "SpectralContent") {
            std::string typeStr;
            fin >> typeStr;

            if (typeStr == "Monochrome") {
                colorModel = ColorModel::MONOCHROME;
                numChannels = 1;
            }
            else if (typeStr == "XYZ") {
                colorModel = ColorModel::XYZ;
                numChannels = 3;
            }
        }
        else if (headStr == "ScatterAzimuth") {
            int numOutPhi;
            fin >> numOutPhi;
            for (int i = 0; i < numOutPhi; ++i) {
                float angle;
                fin >> angle;
                outPhiDegrees.push_back(angle);
            }
        }
        else if (headStr == "ScatterRadial") {
            int numOutTheta;
            fin >> numOutTheta;
            for (int i = 0; i < numOutTheta; ++i) {
                float angle;
                fin >> angle;
                outThetaDegrees.push_back(angle);
            }
        }
        else if (headStr == "DataBegin") {
            break;
        }
    }

    std::vector<Data*> frontBrdfData;
    std::vector<Data*> frontBtdfData;
    std::vector<Data*> backBrdfData;
    std::vector<Data*> backBtdfData;

    int numOutDirSamples = outThetaDegrees.size() * outPhiDegrees.size();

    ignoreCommentLines(fin);

    // Read data.
    Data* data = 0;
    std::string dataStr;
    while (fin >> dataStr) {
        ignoreCommentLines(fin);

        if (!data) {
            data = new Data;
        }

        if (dataStr.empty()) {
            continue;
        }
        else if (dataStr == "AOI") {
            float aoi;
            fin >> aoi;
            if (aoi > 90.0f) {
                aoi = 0;
            }
            data->aoi = aoi;
        }
        else if (dataStr == "POI") {
            float poi;
            fin >> poi;
            data->poi = poi;
        }
        else if (dataStr == "Side") {
            std::string sideStr;
            fin >> sideStr;
            
            if (sideStr == "Front") {
                data->sideType = FRONT;
            }
            else if (sideStr == "Back") {
                data->sideType = BACK;
            }
        }
        else if (dataStr == "Wavelength") {
            float wl;
            fin >> wl;
            data->wavelength = wl;
        }
        else if (dataStr == "ScatterType") {
            std::string scatterStr;
            fin >> scatterStr;
            
            if (scatterStr == "BRDF") {
                data->scatterType = BRDF;
            }
            else if (scatterStr == "BTDF") {
                data->scatterType = BTDF;
            }
        }
        else if (dataStr == "TristimulusValue") {
            std::string trisStr;
            fin >> trisStr;
            
            if (trisStr == "TrisX") {
                data->tristimulusValueType = TRIS_X;
            }
            else if (trisStr == "TrisY") {
                data->tristimulusValueType = TRIS_Y;
            }
            else if (trisStr == "TrisZ") {
                data->tristimulusValueType = TRIS_Z;
            }
        }
        else if (dataStr == "TIS") {
            float tis;
            fin >> tis;
            data->tis = tis;

            ignoreCommentLines(fin);

            data->samples.resize(numOutDirSamples);
            for (int i = 0; i < numOutDirSamples; ++i) {
                float val;
                fin >> val;
                data->samples[i] = val;            
            }

            if (data->sideType == FRONT) {
                if (data->scatterType == BRDF) {
                    frontBrdfData.push_back(data);
                }
                else {
                    frontBtdfData.push_back(data);
                }
            }
            else {
                if (data->scatterType == BRDF) {
                    backBrdfData.push_back(data);
                }
                else {
                    backBtdfData.push_back(data);
                }
            }

            data = 0;
        }
        else if (dataStr == "DataEnd") {
            break;
        }
    }

    SphericalCoordinatesBrdf* frontBrdf = createBrdf(frontBrdfData, outThetaDegrees, outPhiDegrees, colorModel);
    SphericalCoordinatesBrdf* frontBtdf = createBrdf(frontBtdfData, outThetaDegrees, outPhiDegrees, colorModel);
    SphericalCoordinatesBrdf* backBrdf  = createBrdf(backBrdfData,  outThetaDegrees, outPhiDegrees, colorModel);
    SphericalCoordinatesBrdf* backBtdf  = createBrdf(backBtdfData,  outThetaDegrees, outPhiDegrees, colorModel);

    for (auto it = frontBrdfData.begin(); it != frontBrdfData.end(); ++it) { delete *it; }
    for (auto it = frontBtdfData.begin(); it != frontBtdfData.end(); ++it) { delete *it; }
    for (auto it = backBrdfData.begin();  it != backBrdfData.end();  ++it) { delete *it; }
    for (auto it = backBtdfData.begin();  it != backBtdfData.end();  ++it) { delete *it; }

    if (!frontBrdf && !frontBtdf && !backBrdf && !backBtdf) return 0;

    Btdf* fBtdf = frontBtdf ? new Btdf(frontBtdf) : 0;
    Btdf* bBtdf = backBtdf  ? new Btdf(backBtdf)  : 0;
    Bsdf* frontBsdf = new Bsdf(frontBrdf, fBtdf);
    Bsdf* backBsdf  = new Bsdf(backBrdf,  bBtdf);
    Material* frontMaterial = new Material(frontBsdf, 0, 0);
    Material* backMaterial  = new Material(backBsdf,  0, 0);
    TwoSidedMaterial* material = new TwoSidedMaterial(frontMaterial, backMaterial);

    return material;
}

LightToolsBsdfReader::Data::Data() : aoi(0.0f),
                                     poi(0.0f),
                                     wavelength(0.0f),
                                     tis(0.0f),
                                     sideType(FRONT),
                                     scatterType(BRDF) {}

bool LightToolsBsdfReader::Data::cmp(Data* lhs, Data* rhs)
{
    if (lhs->aoi == rhs->aoi) {
        return lhs->tristimulusValueType < rhs->tristimulusValueType;
    }
    else {
        return lhs->aoi < rhs->aoi;
    }
}

SphericalCoordinatesBrdf* LightToolsBsdfReader::createBrdf(std::vector<Data*>&          brdfData,
                                                           const std::vector<float>&    outThetaDegrees,
                                                           const std::vector<float>&    outPhiDegrees,
                                                           ColorModel::Type             colorModel)
{
    if (brdfData.empty()) return 0;

    int numChannels;
    if (colorModel == ColorModel::MONOCHROME) {
        numChannels = 1;
    }
    else if (colorModel == ColorModel::XYZ) {
        numChannels = 3;
    }
    else {
        return 0;
    }

    std::sort(brdfData.begin(), brdfData.end(), Data::cmp);

    std::set<float> inThetaDegrees;
    for (auto it = brdfData.begin(); it != brdfData.end(); ++it) {
        inThetaDegrees.insert((*it)->aoi);
    }

    if (brdfData.size() != inThetaDegrees.size() * numChannels) {
        std::cerr << "[LightToolsBsdfReader::createBrdf] Invalid format." << std::endl;
        return 0;
    }

    SphericalCoordinatesBrdf* brdf = new SphericalCoordinatesBrdf(inThetaDegrees.size(), 1,
                                                                  outThetaDegrees.size(), outPhiDegrees.size(),
                                                                  colorModel,
                                                                  numChannels);
    
    SampleSet* ss = brdf->getSampleSet();

    copy(inThetaDegrees,  ss->getAngles0());
    brdf->setInPhi(0, 0.0f);
    copy(outThetaDegrees, ss->getAngles2());
    copy(outPhiDegrees,   ss->getAngles3());

    ss->getAngles0() = toRadians(ss->getAngles0());
    ss->getAngles2() = toRadians(ss->getAngles2());
    ss->getAngles3() = toRadians(ss->getAngles3());

    for (int i = 0; i < static_cast<int>(brdfData.size()); ++i) {
        Data* data = brdfData.at(i);

        int channelIndex;
        if (colorModel == ColorModel::MONOCHROME) {
            channelIndex = 0;
        }
        else if (colorModel == ColorModel::XYZ) {
            channelIndex = data->tristimulusValueType;
        }
        else {
            return 0;
        }

        ss->setWavelength(channelIndex, data->wavelength);

        int inThIndex = i / numChannels;
        int index = 0;
        for (int outPhIndex = 0; outPhIndex < brdf->getNumOutPhi();   ++outPhIndex)          {
        for (int outThIndex = 0; outThIndex < brdf->getNumOutTheta(); ++outThIndex, ++index) {
            Spectrum& sp = brdf->getSpectrum(inThIndex, 0, outThIndex, outPhIndex);
            sp[channelIndex] = data->samples[index];
        }}

        std::cout << "[LightToolsBsdfReader::read] TIS(inThIndex: " << inThIndex << "): " << data->tis << std::endl;
    }

    // An incoming azimuthal angle of an isotropic LightTools BSDF is 90 degrees.
    SphericalCoordinatesBrdf* rotatedBrdf = rotateOutPhi(*brdf, -PI_2_F);
    rotatedBrdf->clampAngles();

    delete brdf;
    return rotatedBrdf;
}

SphericalCoordinatesBrdf* LightToolsBsdfReader::rotateOutPhi(const SphericalCoordinatesBrdf& brdf, float rotationAngle)
{
    assert(rotationAngle > -2.0f * PI_F && rotationAngle < 2.0f * PI_F);

    if (rotationAngle < 0.0f) {
        rotationAngle += 2.0f * PI_F;
    }

    SphericalCoordinatesBrdf* rotatedBrdf = new SphericalCoordinatesBrdf(brdf);
    SampleSet* ss = rotatedBrdf->getSampleSet();

    ss->checkEqualIntervalAngles();
    if (!ss->isEqualIntervalAngles3()) {
        for (int i = 0; i < rotatedBrdf->getNumOutPhi(); ++i) {
            float outPhi = rotatedBrdf->getOutPhi(i) + rotationAngle;
            if (outPhi > 2.0f * PI_F) {
                outPhi -= 2.0f * PI_F;
            }

            rotatedBrdf->setOutPhi(i, outPhi);
        }

        Arrayf& outPhiAngles = ss->getAngles3();
        std::sort(outPhiAngles.data(), outPhiAngles.data() + outPhiAngles.size());
    }

    for (int inThIndex  = 0; inThIndex  < rotatedBrdf->getNumInTheta();  ++inThIndex)  {
    for (int inPhIndex  = 0; inPhIndex  < rotatedBrdf->getNumInPhi();    ++inPhIndex)  {
    for (int outThIndex = 0; outThIndex < rotatedBrdf->getNumOutTheta(); ++outThIndex) {
    for (int outPhIndex = 0; outPhIndex < rotatedBrdf->getNumOutPhi();   ++outPhIndex) {
        float inTheta  = rotatedBrdf->getInTheta(inThIndex);
        float inPhi    = rotatedBrdf->getInPhi(inPhIndex);
        float outTheta = rotatedBrdf->getOutTheta(outThIndex);
        float outPhi   = rotatedBrdf->getOutPhi(outPhIndex) - rotationAngle;

        if (outPhi < 0.0f) {
            outPhi += 2.0f * PI_F;
        }

        Spectrum sp = brdf.getSpectrum(inTheta, inPhi, outTheta, outPhi);
        rotatedBrdf->setSpectrum(inThIndex, inPhIndex, outThIndex, outPhIndex, sp);
    }}}}

    return rotatedBrdf;
}
