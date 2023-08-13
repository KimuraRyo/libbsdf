// =================================================================== //
// Copyright (C) 2020-2023 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#include <libbsdf/Reader/SsddReader.h>

#include <fstream>

#include <libbsdf/Brdf/HalfDifferenceCoordinatesBrdf.h>
#include <libbsdf/Brdf/Processor.h>
#include <libbsdf/Brdf/SpecularCoordinatesBrdf.h>
#include <libbsdf/Brdf/SphericalCoordinatesBrdf.h>

#include <libbsdf/Reader/ReaderUtility.h>
#include <libbsdf/Reader/SsddUtility.h>

using namespace lb;

Material* SsddReader::read(const std::string& fileName)
{
    // std::ios_base::binary is used to read line endings of CR+LF and LF.
    std::ifstream ifs(fileName.c_str(), std::ios_base::binary);
    if (ifs.fail()) {
        lbError << "[SsddReader::read] Could not open: " << fileName;
        return 0;
    }

    std::ios_base::sync_with_stdio(false);

    Material* material = new Material(nullptr, nullptr, nullptr);

    reader_utility::ignoreCommentLines(ifs, "#");

    FileInfo fileInfo;
    DataInfo* dataInfo = nullptr;
    DataInfo brdfInfo, btdfInfo, srInfo, stInfo;

    std::string lineStr;
    while (!reader_utility::safeGetline(ifs, lineStr).eof()) {
        if (lineStr.empty()) continue;
        reader_utility::ignoreCommentLines(ifs, "#");

        std::stringstream stream(lineStr);
        std::string propStr;
        stream >> propStr;

        // Read the file header
        if (propStr == ssdd::VERSION) {
            stream >> fileInfo.version;
            lbInfo << "[SsddReader::read] SSDD version: " << fileInfo.version;
        }
        else if (propStr == ssdd::SOFTWARE) {
            stream >> fileInfo.software;
            lbInfo << "[SsddReader::read] Software: " << fileInfo.software;
        }
        else if (propStr == ssdd::API) {
            stream >> fileInfo.api;
            lbInfo << "[SsddReader::read] API: " << fileInfo.api;
        }
        else if (propStr == ssdd::DATE) {
            stream >> fileInfo.date;
            lbInfo << "[SsddReader::read] Date of file creation: " << fileInfo.date;
        }

        // Read meta-data.
        else if (propStr == ssdd::DATA_TYPE) {
            std::string typeStr;
            stream >> typeStr;

            if (typeStr == ssdd::DATA_TYPE_BRDF) {
                dataInfo = &brdfInfo;
                dataInfo->dataType = BRDF_DATA;
            }
            else if (typeStr == ssdd::DATA_TYPE_BTDF) {
                dataInfo = &btdfInfo;
                dataInfo->dataType = BTDF_DATA;
            }
            else if (typeStr == ssdd::DATA_TYPE_SPECULAR_REFLECTANCE) {
                dataInfo = &srInfo;
                dataInfo->dataType = SPECULAR_REFLECTANCE_DATA;
            }
            else if (typeStr == ssdd::DATA_TYPE_SPECULAR_TRANSMITTANCE) {
                dataInfo = &stInfo;
                dataInfo->dataType = SPECULAR_TRANSMITTANCE_DATA;
            }
            else {
                lbError << "[SsddReader::readBrdf] Invalid data type: " << typeStr;
                return nullptr;
            }

            lbInfo << "[SsddReader::read] Data type: " << typeStr;
        }

        else if (!dataInfo) {
            lbError << "[SsddReader::readBrdf] Data type not found.";
            return nullptr;
        }

        else if (propStr == ssdd::COLOR_MODEL) {
            std::string typeStr;
            stream >> typeStr;

            if (typeStr == ssdd::COLOR_MODEL_MONOCHROME) {
                dataInfo->colorModel = MONOCHROMATIC_MODEL;
            }
            else if (typeStr == ssdd::COLOR_MODEL_RGB) {
                dataInfo->colorModel = RGB_MODEL;
            }
            else if (typeStr == ssdd::COLOR_MODEL_XYZ) {
                dataInfo->colorModel = XYZ_MODEL;
            }
            else if (typeStr == ssdd::COLOR_MODEL_SPECTRUM) {
                dataInfo->colorModel = SPECTRAL_MODEL;
            }
            else {
                lbError << "[SsddReader::readBrdf] Invalid color model: " << typeStr;
                return nullptr;
            }
        }
        else if (propStr == ssdd::WAVELENGTH_LIST) {
            dataInfo->wavelengths = getList<float>(stream);
        }
        else if (propStr == ssdd::PARAM_TYPE) {
            stream >> dataInfo->paramType;
        }
        else if (propStr == ssdd::REDUCTION_TYPE) {
            std::string typeStr;
            while (stream >> typeStr) {
                int type;

                if (typeStr == ssdd::REDUCTION_TYPE_BILATERAL_SYMMETRY) {
                    type = (asInteger(dataInfo->reductionType) |
                            asInteger(ReductionType::BILATERAL_SYMMETRY));
                }
                else if (typeStr == ssdd::REDUCTION_TYPE_RECIPROCITY) {
                    type = (asInteger(dataInfo->reductionType) |
                            asInteger(ReductionType::RECIPROCITY));
                }
                else {
                    lbError << "[SsddReader::readBrdf] Invalid reduction type: " << typeStr;
                    return nullptr;
                }

                dataInfo->reductionType = static_cast<ReductionType>(type);
            }
        }
        else if (propStr == ssdd::PARAM0_LIST) {
            dataInfo->params0 = getList<double>(stream);
        }
        else if (propStr == ssdd::PARAM1_LIST) {
            dataInfo->params1 = getList<double>(stream);
        }
        else if (propStr == ssdd::PARAM2_LIST) {
            dataInfo->params2 = getList<double>(stream);
        }
        else if (propStr == ssdd::PARAM3_LIST) {
            dataInfo->params3 = getList<double>(stream);
        }
        else if (propStr == ssdd::PARAM4_LIST) {
            dataInfo->params4 = getList<double>(stream);
        }

        // Read optional meta-data.
        else if (propStr == ssdd::SOURCE_TYPE) {
            std::string typeStr;
            stream >> typeStr;

            if (typeStr == ssdd::SOURCE_TYPE_EDITED) {
                dataInfo->sourceType = EDITED_SOURCE;
            }
            else if (typeStr == ssdd::SOURCE_TYPE_GENERATED) {
                dataInfo->sourceType = GENERATED_SOURCE;
            }
            else if (typeStr == ssdd::SOURCE_TYPE_MEASURED) {
                dataInfo->sourceType = MEASURED_SOURCE;
            }
            else {
                lbError << "[SsddReader::readBrdf] Invalid source type: " << typeStr;
                return nullptr;
            }
        }
        else if (propStr == ssdd::DEVICE) {
            stream >> dataInfo->device;
            lbInfo << "[SsddReader::read] Measurement device: " << dataInfo->device;
        }
        else if (propStr == ssdd::CREATION_DATE) {
            stream >> dataInfo->creationDate;
            lbInfo << "[SsddReader::read] Date of data creation: " << dataInfo->creationDate;
        }
        else if (propStr == ssdd::MEASUREMENT_DATE) {
            stream >> dataInfo->measurementDate;
            lbInfo << "[SsddReader::read] Date of measurement: " << dataInfo->measurementDate;
        }

        // Read tabular data.
        else if (propStr == ssdd::DATA) {
            std::string typeStr;
            stream >> typeStr;
            dataInfo->dataMode = typeStr;

            if (dataInfo->params1.empty()) {
                dataInfo->params1.push_back(0.0f);
            }

            if (dataInfo->dataType == BRDF_DATA) {
                std::shared_ptr<Brdf> brdf = readBrdf(ifs, brdfInfo);
                if (!brdf) {
                    delete material;
                    return nullptr;
                }

                std::shared_ptr<Bsdf> bsdf = material->getBsdf();
                if (bsdf) {
                    bsdf->setBrdf(brdf);
                }
                else {
                    bsdf.reset(new Bsdf(brdf, nullptr));
                    material->setBsdf(bsdf);
                }
            }
            else if (dataInfo->dataType == BTDF_DATA) {
                std::shared_ptr<Brdf> brdf = readBrdf(ifs, *dataInfo);
                if (!brdf) {
                    delete material;
                    return nullptr;
                }

                std::shared_ptr<Bsdf> bsdf = material->getBsdf();
                std::shared_ptr<Btdf> btdf = std::make_shared<Btdf>(brdf);
                if (bsdf) {
                    bsdf->setBtdf(btdf);
                }
                else {
                    bsdf.reset(new Bsdf(nullptr, btdf));
                    material->setBsdf(bsdf);
                }
            }
            else if (dataInfo->dataType == SPECULAR_REFLECTANCE_DATA) {
                std::shared_ptr<SampleSet2D> ss2 = readSpecularReflectances(ifs, *dataInfo);
                if (!ss2) {
                    delete material;
                    return nullptr;
                }

                material->setSpecularReflectances(ss2);
            }
            else if (dataInfo->dataType == SPECULAR_TRANSMITTANCE_DATA) {
                std::shared_ptr<SampleSet2D> ss2 = readSpecularReflectances(ifs, *dataInfo);
                if (!ss2) {
                    delete material;
                    return nullptr;
                }

                material->setSpecularTransmittances(ss2);
            }
        }
        else {
            reader_utility::logNotImplementedKeyword(propStr);
        }
    }

    return material;
}

std::shared_ptr<Brdf> SsddReader::readBrdf(std::ifstream& ifs, const DataInfo& dataInfo)
{
    std::shared_ptr<Brdf> brdf;

    if (dataInfo.paramType == ssdd::PARAM_TYPE_HALF_DIFF) {
        brdf = std::make_shared<HalfDifferenceCoordinatesBrdf>(static_cast<int>(dataInfo.params0.size()),
                                                               static_cast<int>(dataInfo.params1.size()),
                                                               static_cast<int>(dataInfo.params2.size()),
                                                               static_cast<int>(dataInfo.params3.size()),
                                                               dataInfo.colorModel,
                                                               static_cast<int>(dataInfo.wavelengths.size()));
    }
    else if (dataInfo.paramType == ssdd::PARAM_TYPE_SPECULAR) {
        brdf = std::make_shared<SpecularCoordinatesBrdf>(static_cast<int>(dataInfo.params0.size()),
                                                         static_cast<int>(dataInfo.params1.size()),
                                                         static_cast<int>(dataInfo.params2.size()),
                                                         static_cast<int>(dataInfo.params3.size()),
                                                         dataInfo.colorModel,
                                                         static_cast<int>(dataInfo.wavelengths.size()));
    }
    else if (dataInfo.paramType == ssdd::PARAM_TYPE_SPHERICAL) {
        brdf = std::make_shared<SphericalCoordinatesBrdf>(static_cast<int>(dataInfo.params0.size()),
                                                          static_cast<int>(dataInfo.params1.size()),
                                                          static_cast<int>(dataInfo.params2.size()),
                                                          static_cast<int>(dataInfo.params3.size()),
                                                          dataInfo.colorModel,
                                                          static_cast<int>(dataInfo.wavelengths.size()));
    }
    else {
        lbError << "[SsddReader::readBrdf] Invalid parameter type: " << dataInfo.paramType;
        return nullptr;
    }

    brdf->setReductionType(dataInfo.reductionType);

    SampleSet* ss = brdf->getSampleSet();

    Arrayd& angles0 = ss->getAngles0();
    Arrayd& angles1 = ss->getAngles1();
    Arrayd& angles2 = ss->getAngles2();
    Arrayd& angles3 = ss->getAngles3();

    array_util::copy(dataInfo.params0, &angles0);
    array_util::copy(dataInfo.params1, &angles1);
    array_util::copy(dataInfo.params2, &angles2);
    array_util::copy(dataInfo.params3, &angles3);

    angles0 = toRadians(angles0);
    angles1 = toRadians(angles1);
    angles2 = toRadians(angles2);
    angles3 = toRadians(angles3);

    // Set offset angles for lb::SpecularCoordinatesBrdf.
    auto specBrdf = dynamic_cast<SpecularCoordinatesBrdf*>(brdf.get());
    if (specBrdf &&
        dataInfo.params4.size() == dataInfo.params0.size()) {
        specBrdf->getSpecularOffsets().resize(dataInfo.params4.size());
        Arrayd& offsets = specBrdf->getSpecularOffsets();
        array_util::copy(dataInfo.params4, &offsets);
        offsets = toRadians(offsets);
    }

    array_util::copy(dataInfo.wavelengths, &ss->getWavelengths());

    if (dataInfo.dataMode == ssdd::DATA_ASCII) {
        if (!readAsciiData(ifs, brdf->getSampleSet())) {
            return nullptr;
        }
    }
    else if (dataInfo.dataMode == ssdd::DATA_BINARY) {
        if (!readBinaryData(ifs, brdf->getSampleSet())) {
            return nullptr;
        }
    }
    else {
        lbError << "[SsddReader::readBrdf] Invalid data mode: " << dataInfo.dataMode;
        return nullptr;
    }

    if (!ss->validate()) {
        lbError << "[SsddReader::readBrdf] Invalid data.";
        return nullptr;
    }

    brdf->clampAngles();

    HalfDifferenceCoordinatesBrdf* hdBrdf = dynamic_cast<HalfDifferenceCoordinatesBrdf*>(brdf.get());
    if (hdBrdf && hasSameEnumerator(brdf->getReductionType(), ReductionType::RECIPROCITY)) {
        brdf.reset(fillAnglesUsingReciprocity(*hdBrdf));
    }

    if (hasSameEnumerator(brdf->getReductionType(), ReductionType::BILATERAL_SYMMETRY)) {
        brdf.reset(fillAnglesUsingBilateralSymmetry(*brdf));
    }

    return brdf;
}

std::shared_ptr<SampleSet2D> SsddReader::readSpecularReflectances(std::ifstream& ifs, const DataInfo& dataInfo)
{
    std::shared_ptr<SampleSet2D> ss2 = std::make_shared<SampleSet2D>(static_cast<int>(dataInfo.params0.size()),
                                                                     static_cast<int>(dataInfo.params1.size()),
                                                                     dataInfo.colorModel,
                                                                     static_cast<int>(dataInfo.wavelengths.size()));

    array_util::copy(dataInfo.params0, &ss2->getThetaArray());
    array_util::copy(dataInfo.params1, &ss2->getPhiArray());

    ss2->getThetaArray()    = toRadians(ss2->getThetaArray());
    ss2->getPhiArray()      = toRadians(ss2->getPhiArray());

    ss2->clampAngles();

    array_util::copy(dataInfo.wavelengths, &ss2->getWavelengths());

    if (dataInfo.dataMode == ssdd::DATA_ASCII) {
        if (!readAsciiData(ifs, ss2.get())) {
            return nullptr;
        }
    }
    else if (dataInfo.dataMode == ssdd::DATA_BINARY) {
        if (!readBinaryData(ifs, ss2.get())) {
            return nullptr;
        }
    }
    else {
        lbError << "[SsddReader::readSpecularReflectances] Invalid data mode: " << dataInfo.dataMode;
        return nullptr;
    }

    if (!ss2->validate()) {
        lbError << "[SsddReader::readSpecularReflectances] Invalid data.";
        return nullptr;
    }

    return ss2;
}

bool SsddReader::readAsciiData(std::ifstream& ifs, SampleSet* ss)
{
    std::string lineStr;
    size_t index = 0;
    while (!reader_utility::safeGetline(ifs, lineStr).eof()) {
        if (lineStr.empty()) continue;
        reader_utility::ignoreCommentLines(ifs, "#");

        std::stringstream stream(lineStr);
        std::vector<float> values = getList<float>(stream);

        if (values.size() != ss->getNumWavelengths()) {
            lbError << "[SsddReader::readAsciiData] Invalid data found: " << lineStr;
            return false;
        }

        array_util::copy(values, &ss->getSpectra().at(index));

        if (index == ss->getSpectra().size() - 1) break;

        ++index;
    }

    if (index != ss->getSpectra().size() - 1) {
        lbError << "[SsddReader::readAsciiData] Invalid data format.";
        return false;
    }

    return true;
}

bool SsddReader::readBinaryData(std::ifstream& ifs, SampleSet* ss)
{
    for (int i3 = 0; i3 < ss->getNumAngles3(); ++i3) {
    for (int i2 = 0; i2 < ss->getNumAngles2(); ++i2) {
    for (int i1 = 0; i1 < ss->getNumAngles1(); ++i1) {
    for (int i0 = 0; i0 < ss->getNumAngles0(); ++i0) {
        Spectrum sp(ss->getNumWavelengths());

        // to do: Perform a byte swap in the big-endian system.
        if (!ifs.read(reinterpret_cast<char*>(sp.data()), sp.size() * sizeof(typename Spectrum::Scalar))) {
            lbError << "[SsddReader::readBinaryData] Failed to read data.";
            return false;
        }

        ss->setSpectrum(i0, i1, i2, i3, sp);
    }}}}

    return true;
}

bool SsddReader::readAsciiData(std::ifstream& ifs, SampleSet2D* ss2)
{
    std::string lineStr;
    size_t index = 0;
    while (!reader_utility::safeGetline(ifs, lineStr).eof()) {
        if (lineStr.empty()) continue;
        reader_utility::ignoreCommentLines(ifs, "#");

        std::stringstream stream(lineStr);
        std::vector<float> values = getList<float>(stream);

        if (values.size() != ss2->getNumWavelengths()) {
            lbError << "[SsddReader::readAsciiData] Invalid data found: " << lineStr;
            return false;
        }

        array_util::copy(values, &ss2->getSpectra().at(index));

        if (index == ss2->getSpectra().size() - 1) break;

        ++index;
    }

    if (index != ss2->getSpectra().size() - 1) {
        lbError << "[SsddReader::readAsciiData] Invalid data format.";
        return false;
    }

    return true;
}

bool SsddReader::readBinaryData(std::ifstream& ifs, SampleSet2D* ss2)
{
    for (int i1 = 0; i1 < ss2->getNumPhi();   ++i1) {
    for (int i0 = 0; i0 < ss2->getNumTheta(); ++i0) {
        Spectrum sp(ss2->getNumWavelengths());

        // to do: Perform a byte swap in the big-endian system.
        if (!ifs.read(reinterpret_cast<char*>(sp.data()), sp.size() * sizeof(typename Spectrum::Scalar))) {
            lbError << "[SsddReader::readBinaryData] Failed to read data.";
            return false;
        }

        ss2->setSpectrum(i0, i1, sp);
    }}

    return true;
}
