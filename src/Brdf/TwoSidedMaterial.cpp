// =================================================================== //
// Copyright (C) 2015 Kimura Ryo                                       //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#include <libbsdf/Brdf/TwoSidedMaterial.h>

using namespace lb;

TwoSidedMaterial::TwoSidedMaterial(Material* frontMaterial,
                                   Material* backMaterial)
                                   : frontMaterial_(frontMaterial),
                                     backMaterial_(backMaterial) {}

TwoSidedMaterial::~TwoSidedMaterial()
{
    delete frontMaterial_;
    delete backMaterial_;
}
