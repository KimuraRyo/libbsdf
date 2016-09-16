// =================================================================== //
// Copyright (C) 2014-2016 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

/*!
 * \file    ReaderUtility.h
 * \brief   The ReaderUtility.h header file includes the utility functions for file readers.
 */

#ifndef LIBBSDF_READER_UTILITY_H
#define LIBBSDF_READER_UTILITY_H

#include <string>
#include <istream>
#include <iostream>
#include <limits>
#include <algorithm>

#include <libbsdf/Common/Global.h>

namespace lb {
namespace reader_utility {

/*! Outputs an unsupported keyword. */
void logNotImplementedKeyword(const std::string& keyword);

/*! Skips a line. */
void ignoreLine(std::istream& stream);

/*! Skips comment lines. */
void ignoreCommentLines(std::istream& stream, const std::string& lineHead);

/*! Converts a string to lower-case. */
std::string toLower(const std::string& str);

/*! Returns true if the string ends with \a suffix. */
bool hasSuffix(const std::string &fileName, const std::string &suffix);

/*! Classifies the type of a file. */
FileType classifyFile(const std::string& fileName);

} // namespace reader_utility

inline void reader_utility::logNotImplementedKeyword(const std::string& keyword)
{
    std::cerr << "Not implemented: " << "\"" << keyword << "\"" << std::endl;
}

inline void reader_utility::ignoreLine(std::istream& stream)
{
    stream.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
}

inline std::string reader_utility::toLower(const std::string& str)
{
    std::string lowerStr(str);
    std::transform(lowerStr.begin(), lowerStr.end(), lowerStr.begin(), ::tolower);
    return lowerStr;
}

} // namespace lb

#endif // LIBBSDF_READER_UTILITY_H
