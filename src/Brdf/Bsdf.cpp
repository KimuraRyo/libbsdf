// =================================================================== //
// Copyright (C) 2015-2020 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#include <libbsdf/Brdf/Bsdf.h>

using namespace lb;

Bsdf::Bsdf(std::shared_ptr<Brdf> brdf,
           std::shared_ptr<Btdf> btdf)
           : brdf_(brdf),
             btdf_(btdf)
{
    lbTrace << "[Bsdf::Bsdf]";
}

Bsdf::~Bsdf()
{
    lbTrace << "[Bsdf::~Bsdf]";
}

bool Bsdf::validate(bool verbose) const
{
    bool brdfInvalid = (brdf_ && !brdf_->validate(verbose));
    bool btdfInvalid = (btdf_ && !btdf_->validate(verbose));

    return (!brdfInvalid && !btdfInvalid);
}
