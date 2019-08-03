// =================================================================== //
// Copyright (C) 2018-2019 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#include <ArgumentParser.h>

#include <algorithm>

ArgumentParser::ArgumentParser(int argc, char** argv)
{
    for (int i = 1; i < argc; ++i) {
        tokens_.push_back(std::string(argv[i]));
    }
}

bool ArgumentParser::read(const std::string& str)
{
    auto it = std::find(tokens_.begin(), tokens_.end(), str);
    if (it != tokens_.end()) {
        tokens_.erase(it);
        return true;
    }

    return false;
}

bool ArgumentParser::read(const std::string& option, std::string* value)
{
    auto optionIt = std::find(tokens_.begin(), tokens_.end(), option);
    if (optionIt != tokens_.end()) {
        auto valueIt = ++optionIt;
        if (valueIt != tokens_.end()) {
            *value = *valueIt;

            tokens_.erase(valueIt);
            tokens_.erase(std::find(tokens_.begin(), tokens_.end(), option));

            return true;
        }
    }

    return false;
}
