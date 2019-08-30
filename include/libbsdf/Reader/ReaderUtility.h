// =================================================================== //
// Copyright (C) 2014-2019 Kimura Ryo                                  //
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
#include <limits>
#include <algorithm>

#include <libbsdf/Common/Global.h>
#include <libbsdf/Common/Log.h>

namespace lb {
namespace reader_utility {

/*! \brief Outputs an unsupported keyword. */
void logNotImplementedKeyword(const std::string& keyword);

/*! \brief Skips a line. */
void ignoreLine(std::istream& stream);

/*! \brief Skips comment lines. */
void ignoreCommentLines(std::istream& stream, const std::string& lineHead);

/*!
 * \brief \a getline function to handle all three line endings ("\r", "\n" and "\r\n").
 * See https://stackoverflow.com/questions/6089231/getting-std-ifstream-to-handle-lf-cr-and-crlf
*/
std::istream& safeGetline(std::istream& stream, std::string& token);

/*! \brief Converts a string to lower-case. */
std::string toLower(const std::string& str);

/*! \brief Returns true if two strings converted to lower-case are equal. */
bool isEqual(const std::string& lhs, const std::string& rhs);

/*! \brief Returns true if the string ends with \a suffix. */
bool hasSuffix(const std::string &fileName, const std::string &suffix);

/*! \brief Classifies the type of a file. */
FileType classifyFile(const std::string& fileName);

} // namespace reader_utility

/*
 * Implementation
 */

inline void reader_utility::logNotImplementedKeyword(const std::string& keyword)
{
    lbError << "Not implemented: " << "\"" << keyword << "\"";
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

inline bool reader_utility::isEqual(const std::string& lhs, const std::string& rhs)
{
    return (toLower(lhs) == toLower(rhs));
}

} // namespace lb

#endif // LIBBSDF_READER_UTILITY_H
