// =================================================================== //
// Copyright (C) 2014-2015 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

/*!
 * \file    Processor.h
 * \brief   The Processor.h header file includes the functions to process BRDFs.
 */

#ifndef LIBBSDF_PROCESSOR_H
#define LIBBSDF_PROCESSOR_H

#include <libbsdf/Brdf/Brdf.h>

namespace lb {

/*! \brief Divides a BRDF by the cosine of the outgoing polar angle.
 *
 * This function converts from a CCBRDF to a BRDF.
 */
void divideByCosineOutTheta(Brdf* brdf);

} // namespace lb

#endif // LIBBSDF_PROCESSOR_H
