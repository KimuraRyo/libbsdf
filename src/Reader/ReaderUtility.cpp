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
        bool isComment = (static_cast<int>(peekStr.size()) >= strSize &&
                          peekStr.substr(0, strSize) == lineHead);
        if (isComment) {
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
