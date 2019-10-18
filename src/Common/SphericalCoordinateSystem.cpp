// =================================================================== //
// Copyright (C) 2014-2019 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#include <libbsdf/Common/SphericalCoordinateSystem.h>

using namespace lb;

const char SphericalCoordinateSystem::ANGLE0_NAME[] = "Incoming polar angle";
const char SphericalCoordinateSystem::ANGLE1_NAME[] = "Incoming azimuthal angle";
const char SphericalCoordinateSystem::ANGLE2_NAME[] = "Outgoing polar angle";
const char SphericalCoordinateSystem::ANGLE3_NAME[] = "Outgoing azimuthal angle";

const float SphericalCoordinateSystem::MIN_ANGLE0 = 0.0f;
const float SphericalCoordinateSystem::MIN_ANGLE1 = 0.0f;
const float SphericalCoordinateSystem::MIN_ANGLE2 = 0.0f;
const float SphericalCoordinateSystem::MIN_ANGLE3 = 0.0f;

const float SphericalCoordinateSystem::MAX_ANGLE0 = decrease(PI_2_F);
const float SphericalCoordinateSystem::MAX_ANGLE1 = decrease(TAU_F);
const float SphericalCoordinateSystem::MAX_ANGLE2 = decrease(PI_2_F);
const float SphericalCoordinateSystem::MAX_ANGLE3 = decrease(TAU_F);
