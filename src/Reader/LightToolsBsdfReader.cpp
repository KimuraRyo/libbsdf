// =================================================================== //
// Copyright (C) 2015-2019 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#include <libbsdf/Reader/LightToolsBsdfReader.h>

#include <fstream>
#include <set>

#include <libbsdf/Brdf/Processor.h>

using namespace lb;

TwoSidedMaterial* LightToolsBsdfReader::read(const std::string& fileName)
{
    // std::ios_base::binary is used to read line endings of CR+LF and LF.
    std::ifstream ifs(fileName.c_str(), std::ios_base::binary);
    if (ifs.fail()) {
        lbError << "[LightToolsBsdfReader::read] Could not open: " << fileName;
        return 0;
    }

    std::ios_base::sync_with_stdio(false);

    SymmetryType symmetryType = UNKNOWN_SYMMETRY;
    ColorModel colorModel = UNKNOWN_MODEL;

    std::vector<float> outThetaDegrees;
    std::vector<float> outPhiDegrees;

    int numChannels = 1;

    ignoreCommentLines(ifs);

    // Read a header.
    std::string headStr;
    while (ifs >> headStr) {
        ignoreCommentLines(ifs);

        if (headStr.empty()) {
            continue;
        }
        else if (headStr == "Symmetry") {
            std::string typeStr;
            ifs >> typeStr;

            if (typeStr == "PlaneSymmetrical") {
                symmetryType = PLANE_SYMMETRICAL;
            }
            else if (typeStr == "Asymmetric") {
                symmetryType = ASYMMETRICAL;
            }
            else {
                reader_utility::logNotImplementedKeyword(typeStr);
                return 0;
            }
        }
        else if (headStr == "SpectralContent") {
            std::string typeStr;
            ifs >> typeStr;

            if (typeStr == "Monochrome") {
                colorModel = MONOCHROMATIC_MODEL;
                numChannels = 1;
            }
            else if (typeStr == "XYZ") {
                colorModel = XYZ_MODEL;
                numChannels = 3;
            }
            else {
                reader_utility::logNotImplementedKeyword(typeStr);
                return 0;
            }
        }
        else if (headStr == "ScatterAzimuth") {
            int numOutPhi;
            ifs >> numOutPhi;
            for (int i = 0; i < numOutPhi; ++i) {
                float angle;
                ifs >> angle;
                outPhiDegrees.push_back(angle);
            }
        }
        else if (headStr == "ScatterRadial") {
            int numOutTheta;
            ifs >> numOutTheta;
            for (int i = 0; i < numOutTheta; ++i) {
                float angle;
                ifs >> angle;
                outThetaDegrees.push_back(angle);
            }
        }
        else if (headStr == "DataBegin") {
            break;
        }
    }

    std::vector<DataBlock*> frontBrdfData;
    std::vector<DataBlock*> frontBtdfData;
    std::vector<DataBlock*> backBrdfData;
    std::vector<DataBlock*> backBtdfData;

    if (outThetaDegrees.empty() ||
        outPhiDegrees.empty()) {
        lbError << "[LightToolsBsdfReader::read] Invalid format.";
        return 0;
    }

    size_t numOutDirSamples = outThetaDegrees.size() * outPhiDegrees.size();

    ignoreCommentLines(ifs);

    // Read data.
    DataBlock* data = 0;
    std::string dataStr;
    while (ifs >> dataStr) {
        ignoreCommentLines(ifs);

        if (!data) {
            data = new DataBlock;
        }

        if (dataStr.empty()) {
            continue;
        }
        else if (dataStr == "AOI") {
            float aoi;
            ifs >> aoi;
            if (aoi > 90.0f) {
                aoi = 0;
            }
            data->aoi = aoi;
        }
        else if (dataStr == "POI") {
            ifs >> data->poi;
        }
        else if (dataStr == "Side") {
            std::string sideStr;
            ifs >> sideStr;

            if (sideStr == "Front") {
                data->sideType = FRONT_SIDE;
            }
            else if (sideStr == "Back") {
                data->sideType = BACK_SIDE;
            }
            else {
                reader_utility::logNotImplementedKeyword(sideStr);
                return 0;
            }
        }
        else if (dataStr == "Wavelength") {
            ifs >> data->wavelength;
        }
        else if (dataStr == "ScatterType") {
            std::string scatterStr;
            ifs >> scatterStr;

            if (scatterStr == "BRDF") {
                data->dataType = BRDF_DATA;
            }
            else if (scatterStr == "BTDF") {
                data->dataType = BTDF_DATA;
            }
            else {
                reader_utility::logNotImplementedKeyword(scatterStr);
                return 0;
            }
        }
        else if (dataStr == "TristimulusValue") {
            std::string trisStr;
            ifs >> trisStr;

            if (trisStr == "TrisX") {
                data->tristimulusValueType = TRIS_X;
            }
            else if (trisStr == "TrisY") {
                data->tristimulusValueType = TRIS_Y;
            }
            else if (trisStr == "TrisZ") {
                data->tristimulusValueType = TRIS_Z;
            }
            else {
                reader_utility::logNotImplementedKeyword(trisStr);
                return 0;
            }
        }
        else if (dataStr == "TIS") {
            ifs >> data->tis;

            ignoreCommentLines(ifs);

            data->samples.resize(numOutDirSamples);
            for (size_t i = 0; i < numOutDirSamples; ++i) {
                std::string valueStr;
                ifs >> valueStr;
                data->samples[i] = static_cast<Arrayf::Scalar>(std::atof(valueStr.c_str()));
            }

            if (data->sideType == FRONT_SIDE) {
                if (data->dataType == BRDF_DATA) {
                    frontBrdfData.push_back(data);
                }
                else {
                    frontBtdfData.push_back(data);
                }
            }
            else {
                if (data->dataType == BRDF_DATA) {
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

    std::shared_ptr<Brdf> frontBrdf(createBrdf(frontBrdfData, outThetaDegrees, outPhiDegrees, colorModel));
    std::shared_ptr<Brdf> frontBtdf(createBrdf(frontBtdfData, outThetaDegrees, outPhiDegrees, colorModel));
    std::shared_ptr<Brdf> backBrdf (createBrdf(backBrdfData,  outThetaDegrees, outPhiDegrees, colorModel));
    std::shared_ptr<Brdf> backBtdf (createBrdf(backBtdfData,  outThetaDegrees, outPhiDegrees, colorModel));

    for (auto it = frontBrdfData.begin(); it != frontBrdfData.end(); ++it) { delete *it; }
    for (auto it = frontBtdfData.begin(); it != frontBtdfData.end(); ++it) { delete *it; }
    for (auto it = backBrdfData.begin();  it != backBrdfData.end();  ++it) { delete *it; }
    for (auto it = backBtdfData.begin();  it != backBtdfData.end();  ++it) { delete *it; }

    if (!frontBrdf && !frontBtdf && !backBrdf && !backBtdf) return 0;

    std::shared_ptr<Btdf> fBtdf = frontBtdf ? std::make_shared<Btdf>(frontBtdf) : 0;
    std::shared_ptr<Btdf> bBtdf = backBtdf  ? std::make_shared<Btdf>(backBtdf)  : 0;

    std::shared_ptr<Bsdf> frontBsdf = std::make_shared<Bsdf>(frontBrdf, fBtdf);
    std::shared_ptr<Bsdf> backBsdf  = std::make_shared<Bsdf>(backBrdf,  bBtdf);

    std::shared_ptr<Material> frontMaterial = std::make_shared<Material>(frontBsdf, nullptr, nullptr);
    std::shared_ptr<Material> backMaterial  = std::make_shared<Material>(backBsdf, nullptr, nullptr);

    TwoSidedMaterial* material = new TwoSidedMaterial(frontMaterial, backMaterial);

    return material;
}

LightToolsBsdfReader::DataBlock::DataBlock() : aoi(0.0f),
                                               poi(0.0f),
                                               wavelength(0.0f),
                                               tis(0.0f),
                                               sideType(FRONT_SIDE),
                                               dataType(BRDF_DATA) {}

bool LightToolsBsdfReader::DataBlock::cmp(DataBlock* lhs, DataBlock* rhs)
{
    if (lhs->aoi == rhs->aoi) {
        return lhs->tristimulusValueType < rhs->tristimulusValueType;
    }
    else {
        return lhs->aoi < rhs->aoi;
    }
}

SphericalCoordinatesBrdf* LightToolsBsdfReader::createBrdf(std::vector<DataBlock*>&     brdfData,
                                                           const std::vector<float>&    outThetaDegrees,
                                                           const std::vector<float>&    outPhiDegrees,
                                                           ColorModel                   colorModel)
{
    if (brdfData.empty()) return 0;

    int numChannels;
    if (colorModel == MONOCHROMATIC_MODEL) {
        numChannels = 1;
    }
    else if (colorModel == XYZ_MODEL) {
        numChannels = 3;
    }
    else {
        return 0;
    }

    std::sort(brdfData.begin(), brdfData.end(), DataBlock::cmp);

    std::set<float> inThetaDegrees;
    for (auto it = brdfData.begin(); it != brdfData.end(); ++it) {
        inThetaDegrees.insert((*it)->aoi);
    }

    if (brdfData.size() != inThetaDegrees.size() * numChannels) {
        lbError << "[LightToolsBsdfReader::createBrdf] Invalid format.";
        return 0;
    }

    SphericalCoordinatesBrdf* brdf = new SphericalCoordinatesBrdf(static_cast<int>(inThetaDegrees.size()),
                                                                  1,
                                                                  static_cast<int>(outThetaDegrees.size()),
                                                                  static_cast<int>(outPhiDegrees.size()),
                                                                  colorModel,
                                                                  numChannels);
    SampleSet* ss = brdf->getSampleSet();

    copyArray(inThetaDegrees,  &ss->getAngles0());
    brdf->setInPhi(0, 0.0f);
    copyArray(outThetaDegrees, &ss->getAngles2());
    copyArray(outPhiDegrees,   &ss->getAngles3());

    ss->getAngles0() = toRadians(ss->getAngles0());
    ss->getAngles2() = toRadians(ss->getAngles2());
    ss->getAngles3() = toRadians(ss->getAngles3());

    for (int i = 0; i < static_cast<int>(brdfData.size()); ++i) {
        DataBlock* data = brdfData.at(i);

        int channelIndex;
        if (colorModel == MONOCHROMATIC_MODEL) {
            channelIndex = 0;
        }
        else if (colorModel == XYZ_MODEL) {
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

        lbInfo << "[LightToolsBsdfReader::read] TIS(inThIndex: " << inThIndex << "): " << data->tis;
    }

    // An incoming azimuthal angle of an isotropic LightTools BSDF is 90 degrees.
    SphericalCoordinatesBrdf* rotatedBrdf = rotateOutPhi(*brdf, -PI_2_F);
    rotatedBrdf->clampAngles();
    brdf->setSourceType(MEASURED_SOURCE);

    SampleSet* rotSs = rotatedBrdf->getSampleSet();
    rotSs->updateAngleAttributes();
    if (rotSs->isOneSide()) {
        SphericalCoordinatesBrdf* filledBrdf = fillSymmetricBrdf(rotatedBrdf);
        filledBrdf->expandAngles(false, false, false, true);

        delete rotatedBrdf;
        rotatedBrdf = filledBrdf;
    }

    delete brdf;
    return rotatedBrdf;
}
