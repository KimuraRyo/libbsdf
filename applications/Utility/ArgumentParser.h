// =================================================================== //
// Copyright (C) 2018-2019 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#ifndef LIBBSDF_APPLICATIONS_ARGUMENT_PARSER_H
#define LIBBSDF_APPLICATIONS_ARGUMENT_PARSER_H

#include <iostream>
#include <string>
#include <vector>

/*!
 * \class   ArgumentParser
 * \brief   The ArgumentParser class provides functions to parse command-line arguments in C++.
 */
class ArgumentParser
{
public:
    /*! Results of parameter acquisition. */
    enum ResultType {
        OK,
        NOT_FOUND,
        ERROR
    };

    /*! Constructs tokens from command-line arguments. */
    ArgumentParser(int argc, char** argv);

    /*!
     * Returns true, if \a str is found.
     * \a str is removed from tokens.
     */
    bool read(const std::string& str);

    /*!
     * Returns true, if \a option and \a value are found.
     * \a option and \a value are removed from tokens.
     */
    bool read(const std::string& option, std::string* value);

    /*!
     * Returns true, if \a option and \a value are found.
     * \a option and \a value are removed from tokens.
     */
    template <typename T>
    ResultType read(const std::string& option, T* value);

    /*! Gets the list of tokens. */
    std::vector<std::string>& getTokens() { return tokens_; }

    /*! Validates the number of tokens. If it is invalid, tokens are outputted to std::cerr. */
    bool validate(int numTokens);

private:
    std::vector<std::string> tokens_;
};

template <typename T>
ArgumentParser::ResultType ArgumentParser::read(const std::string& option, T* value)
{
    std::string valStr;
    if (read(option, &valStr)) {
        char* end;
        *value = static_cast<T>(std::strtod(valStr.c_str(), &end));
        if (*end != '\0') {
            std::cerr << "Invalid value (" << option << "): " << valStr << std::endl;
            return ERROR;
        }

        return OK;
    }

    return NOT_FOUND;
}

#endif // LIBBSDF_APPLICATIONS_ARGUMENT_PARSER_H
