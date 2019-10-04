// =================================================================== //
// Copyright (C) 2019 Kimura Ryo                                       //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#include <Utility.h>

#include <libbsdf/Common/Version.h>

using namespace lb;

void app_utility::showAppVersion(const std::string& name, const std::string& version)
{
    std::cout
        << "Version: " << name << " " << version
        << " (libbsdf-" << getVersion() << ")" << std::endl;
}

std::string app_utility::createComments(int                 argc,
                                        char**              argv,
                                        const std::string&  name,
                                        const std::string&  version)
{
    std::string comments("Software: " + name + "-" + version);
    comments += "\n;; Arguments:";
    for (int i = 1; i < argc; ++i) {
        comments += " " + std::string(argv[i]);
    }

    return comments;
}
