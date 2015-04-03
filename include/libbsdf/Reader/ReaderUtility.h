// =================================================================== //
// Copyright (C) 2014-2015 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

/*!
 * \file    ReaderUtility.h
 * \brief   Utility functions for file readers.
 */

#ifndef LIBBSDF_READER_UTILITY_H
#define LIBBSDF_READER_UTILITY_H

#include <string>
#include <fstream>
#include <iostream>
#include <limits>
#include <algorithm>

namespace lb {
namespace reader_utility {

/*! Outputs an unsupported keyword. */
void logNotImplementedKeyword(const std::string& keyword);

/*! Skips a line. */
void ignoreLine(std::ifstream& fin);

/*! Skips comment lines. */
void ignoreCommentLines(std::ifstream& fin, const std::string& lineHead);

/*! Converts a string to lower-case. */
std::string toLower(const std::string& str);

} // namespace reader_utility

inline void reader_utility::logNotImplementedKeyword(const std::string& keyword)
{
    std::cerr << "Not implemented: " << "\"" << keyword << "\"" << std::endl;
}

inline void reader_utility::ignoreLine(std::ifstream& fin)
{
    fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
}

inline std::string reader_utility::toLower(const std::string& str)
{
    std::string lowerStr(str);
    std::transform(lowerStr.begin(), lowerStr.end(), lowerStr.begin(), ::tolower);
    return lowerStr;
}

} // namespace lb

#endif // LIBBSDF_READER_UTILITY_H
