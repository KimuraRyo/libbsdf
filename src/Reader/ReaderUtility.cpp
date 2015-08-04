// =================================================================== //
// Copyright (C) 2014-2015 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#include <libbsdf/Reader/ReaderUtility.h>

using namespace lb;

void reader_utility::ignoreCommentLines(std::ifstream& fin, const std::string& lineHead)
{
    if (fin.fail() || fin.eof()) return;

    std::ifstream::pos_type pos = fin.tellg();

    std::string peekStr;
    while (fin >> peekStr) {
        int strSize = lineHead.size();
        bool comment = (static_cast<int>(peekStr.size()) >= strSize &&
                        peekStr.substr(0, strSize) == lineHead);
        if (comment) {
            ignoreLine(fin);
            pos = fin.tellg();
        }
        else {
            fin.seekg(pos, std::ios_base::beg);
            return;
        }
    }

    if (fin.eof()) {
        fin.clear();
        fin.seekg(pos, std::ios_base::beg);
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
    std::ifstream fin(fileName.c_str());
    if (fin.fail()) {
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
        ignoreCommentLines(fin, "#");

        // Distinguish between LightTools and Zemax.
        std::string token;
        while (fin >> token) {
            ignoreCommentLines(fin, "#");

            if (token == "AOI") {
                return LIGHTTOOLS_FILE;
            }
            else if (token == "AngleOfIncidence") {
                return ZEMAX_FILE;
            }
            else {
                ignoreLine(fin);
            }
        }
    }
    else if (hasSuffix(fileName, ".binary")) {
        return MERL_BINARY_FILE;
    }

    return UNKNOWN_FILE;
}
