// =================================================================== //
// Copyright (C) 2015-2019 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

/*!
 * \file    ReflectanceModelUtility.h
 * \brief   The ReflectanceModelUtility.h header file includes the functions for reflectance models.
 */

#ifndef LIBBSDF_REFLECTANCE_MODEL_UTILITY_H
#define LIBBSDF_REFLECTANCE_MODEL_UTILITY_H

#include <libbsdf/Brdf/Brdf.h>
#include <libbsdf/ReflectanceModel/ReflectanceModel.h>

namespace lb {
namespace reflectance_model_utility {

/*! Sets up a lb::Brdf using an analytic reflectance or transmittance model. */
bool setupTabularBrdf(const ReflectanceModel&   model,
                      Brdf*                     brdf,
                      DataType                  dataType = BRDF_DATA,
                      float                     maxValue = 10000000000.0f);

} // namespace lb
} // namespace reflectance_model_utility

#endif // LIBBSDF_REFLECTANCE_MODEL_UTILITY_H
