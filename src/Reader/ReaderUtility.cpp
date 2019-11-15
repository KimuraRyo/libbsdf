// =================================================================== //
// Copyright (C) 2014-2019 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#include <libbsdf/Reader/ReaderUtility.h>

#include <fstream>

#include <libbsdf/Reader/AstmReader.h>
#include <libbsdf/Reader/DdrReader.h>
#include <libbsdf/Reader/LightToolsBsdfReader.h>
#include <libbsdf/Reader/MerlBinaryReader.h>
#include <libbsdf/Reader/ZemaxBsdfReader.h>

using namespace lb;

void reader_utility::ignoreCommentLines(std::istream& stream, const std::string& lineHead)
{
    if (stream.fail() || stream.eof()) return;

    std::istream::pos_type pos = stream.tellg();

    std::string peekStr;
    while (stream >> peekStr) {
        size_t strSize = lineHead.size();
        bool commentFound = (static_cast<int>(peekStr.size()) >= strSize &&
                             peekStr.substr(0, strSize) == lineHead);
        if (commentFound) {
            ignoreLine(stream);
            pos = stream.tellg();
        }
        else {
            stream.seekg(pos, std::ios_base::beg);
            return;
        }
    }

    if (stream.eof()) {
        stream.clear();
        stream.seekg(pos, std::ios_base::beg);
    }
}

std::istream& reader_utility::safeGetline(std::istream& stream, std::string& token)
{
    token.clear();

    std::istream::sentry se(stream, true);
    std::streambuf* sb = stream.rdbuf();

    for (;;) {
        int c = sb->sbumpc();
        switch (c) {
            case '\n':
                return stream;
            case '\r':
                if (sb->sgetc() == '\n')
                    sb->sbumpc();
                return stream;
            case EOF:
                if (token.empty())
                    stream.setstate(std::ios::eofbit);
                return stream;
            default:
                token += (char)c;
        }
    }
}

bool reader_utility::hasSuffix(const std::string &fileName, const std::string &suffix)
{
    if (fileName.size() >= suffix.size()) {
        return (fileName.compare(fileName.size() - suffix.size(), suffix.size(), suffix) == 0);
    }
    else {
        return false;
    }
}

FileType reader_utility::classifyFile(const std::string& fileName)
{
    // std::ios_base::binary is used to read line endings of CR+LF and LF.
    std::ifstream ifs(fileName.c_str(), std::ios_base::binary);
    if (ifs.fail()) {
        lbError << "[reader_utility::classifyFile] Could not open: " << fileName;
        return UNKNOWN_FILE;
    }

    std::string name = toLower(fileName);

    if (hasSuffix(name, ".astm")) {
        return ASTM_FILE;
    }
    else if (hasSuffix(name, ".csv")) {
        return GCMS4_FILE;
    }
    else if (hasSuffix(name, ".ddr")) {
        return INTEGRA_DDR_FILE;
    }
    else if (hasSuffix(name, ".ddt")) {
        return INTEGRA_DDT_FILE;
    }
    else if (hasSuffix(name, ".sdr")) {
        return INTEGRA_SDR_FILE;
    }
    else if (hasSuffix(name, ".sdt")) {
        return INTEGRA_SDT_FILE;
    }
    else if (hasSuffix(name, ".bsdf")) {
        ignoreCommentLines(ifs, "#");

        // Distinguish between LightTools and Zemax.
        std::string token;
        while (ifs >> token) {
            ignoreCommentLines(ifs, "#");

            if (token == "AOI") {
                return LIGHTTOOLS_FILE;
            }
            else if (token == "AngleOfIncidence") {
                return ZEMAX_FILE;
            }
            else {
                ignoreLine(ifs);
            }
        }
    }
    else if (hasSuffix(name, ".binary")) {
        return MERL_BINARY_FILE;
    }

    return UNKNOWN_FILE;
}

std::shared_ptr<Brdf> reader_utility::read(const std::string&   fileName,
                                           FileType*            fileType,
                                           DataType*            dataType)
{
    std::ifstream ifs(fileName.c_str());
    if (ifs.fail()) {
        lbError << "[reader_utility::read] Could not open: " << fileName;
        return nullptr;
    }

    *fileType = reader_utility::classifyFile(fileName);
    *dataType = UNKNOWN_DATA;

    // Load a BRDF/BTDF.
    std::shared_ptr<Brdf> brdf;
    switch (*fileType) {
        case ASTM_FILE:
            brdf.reset(AstmReader::read(fileName));
            break;
        case INTEGRA_DDR_FILE:
            brdf.reset(DdrReader::read(fileName));
            *dataType = BRDF_DATA;
            break;
        case INTEGRA_DDT_FILE:
            brdf.reset(DdrReader::read(fileName));
            *dataType = BTDF_DATA;
            break;
        case lb::LIGHTTOOLS_FILE: {
            std::unique_ptr<TwoSidedMaterial> material;
            material.reset(LightToolsBsdfReader::read(fileName));
            std::shared_ptr<Brdf> fBrdf = material->getFrontMaterial()->getBsdf()->getBrdf();
            std::shared_ptr<Btdf> fBtdf = material->getFrontMaterial()->getBsdf()->getBtdf();
            std::shared_ptr<Brdf> bBrdf = material->getBackMaterial()->getBsdf()->getBrdf();
            std::shared_ptr<Btdf> bBtdf = material->getBackMaterial()->getBsdf()->getBtdf();

            if      (fBrdf) { brdf = fBrdf;            *dataType = BRDF_DATA; }
            else if (fBtdf) { brdf = fBtdf->getBrdf(); *dataType = BTDF_DATA; }
            else if (bBrdf) { brdf = bBrdf;            *dataType = BRDF_DATA; }
            else if (bBtdf) { brdf = bBtdf->getBrdf(); *dataType = BTDF_DATA; }
            break;
        }
        case lb::MERL_BINARY_FILE:
            brdf.reset(MerlBinaryReader::read(fileName));
            *dataType = BRDF_DATA;
            break;
        case lb::ZEMAX_FILE:
            brdf.reset(ZemaxBsdfReader::read(fileName, dataType));
            break;
        default:
            lbError << "[reader_utility::read] Unsupported file type: " << fileType;
            return nullptr;
    }

    return brdf;
}
