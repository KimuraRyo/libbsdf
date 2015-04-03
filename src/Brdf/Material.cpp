// =================================================================== //
// Copyright (C) 2015 Kimura Ryo                                       //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#include <libbsdf/Brdf/Material.h>

using namespace lb;

Material::Material(Bsdf*        bsdf,
                   SampleSet2D* specularReflectances,
                   SampleSet2D* specularTransmittances,
                   SampleSet2D* reflectionTis,
                   SampleSet2D* transmissionTis)
                   : bsdf_(bsdf),
                     specularReflectances_(specularReflectances),
                     specularTransmittances_(specularTransmittances),
                     reflectionTis_(reflectionTis),
                     transmissionTis_(transmissionTis) {}

Material::~Material()
{
    delete bsdf_;
    delete specularReflectances_;
    delete specularTransmittances_;
    delete reflectionTis_;
    delete transmissionTis_;
}
