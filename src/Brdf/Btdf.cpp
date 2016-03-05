// =================================================================== //
// Copyright (C) 2014-2015 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#include <libbsdf/Brdf/Btdf.h>

using namespace lb;

Btdf::Btdf(Brdf* brdf) : brdf_(brdf) {}

Btdf::~Btdf()
{
    delete brdf_;
}
