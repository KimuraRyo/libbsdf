// =================================================================== //
// Copyright (C) 2015-2020 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#include <libbsdf/Brdf/Material.h>

using namespace lb;

Material::Material(std::shared_ptr<Bsdf>        bsdf,
                   std::shared_ptr<SampleSet2D> specularReflectances,
                   std::shared_ptr<SampleSet2D> specularTransmittances)
                   : bsdf_(bsdf),
                     specularReflectances_(specularReflectances),
                     specularTransmittances_(specularTransmittances) {}

Material::~Material() {}

bool Material::validate(bool verbose) const
{
    bool bsdfInvalid    = (bsdf_                    && !bsdf_->validate(verbose));
    bool srInvalid      = (specularReflectances_    && !specularReflectances_->validate(verbose));
    bool stInvalid      = (specularTransmittances_  && !specularTransmittances_->validate(verbose));

    return (!bsdfInvalid && !srInvalid && !stInvalid);
}
