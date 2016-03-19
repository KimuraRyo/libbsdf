// =================================================================== //
// Copyright (C) 2014-2016 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#include <libbsdf/Reader/ReaderUtility.h>

#include <fstream>

using namespace lb;

void reader_utility::ignoreCommentLines(std::istream& stream, const std::string& lineHead)
{
    if (stream.fail() || stream.eof()) return;

    std::istream::pos_type pos = stream.tellg();

    std::string peekStr;
    while (stream >> peekStr) {
        int strSize = lineHead.size();
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

bool reader_utility::hasSuffix(const std::string &str, const std::string &suffix)
{
    if (str.size() >= suffix.size()) {
        return (str.compare(str.size() - suffix.size(), suffix.size(), suffix) == 0);
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
        std::cerr << "[reader_utility::classifyFile] Could not open: " << fileName << std::endl;
        return UNKNOWN_FILE;
    }

    if (hasSuffix(fileName, ".astm")) {
        return ASTM_FILE;
    }
    else if (hasSuffix(fileName, ".ddr")) {
        return INTEGRA_DDR_FILE;
    }
    else if (hasSuffix(fileName, ".ddt")) {
        return INTEGRA_DDT_FILE;
    }
    else if (hasSuffix(fileName, ".sdr")) {
        return INTEGRA_SDR_FILE;
    }
    else if (hasSuffix(fileName, ".sdt")) {
        return INTEGRA_SDT_FILE;
    }
    else if (hasSuffix(fileName, ".bsdf")) {
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
    else if (hasSuffix(fileName, ".binary")) {
        return MERL_BINARY_FILE;
    }

    return UNKNOWN_FILE;
}
