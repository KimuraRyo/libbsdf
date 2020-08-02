// =================================================================== //
// Copyright (C) 2015-2020 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#include <libbsdf/Brdf/TwoSidedMaterial.h>

using namespace lb;

TwoSidedMaterial::TwoSidedMaterial(std::shared_ptr<Material> frontMaterial,
                                   std::shared_ptr<Material> backMaterial)
                                   : frontMaterial_(frontMaterial),
                                     backMaterial_(backMaterial) {}

TwoSidedMaterial::~TwoSidedMaterial() {}

bool TwoSidedMaterial::validate(bool verbose) const
{
    bool frontMatInvalid    = (frontMaterial_   && !frontMaterial_->validate(verbose));
    bool backMatInvalid     = (backMaterial_    && !backMaterial_->validate(verbose));

    return (!frontMatInvalid && !backMatInvalid);
}
