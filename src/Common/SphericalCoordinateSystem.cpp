// =================================================================== //
// Copyright (C) 2014-2015 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#include <libbsdf/Common/SphericalCoordinateSystem.h>

using namespace lb;

const std::string SphericalCoordinateSystem::ANGLE0_NAME = "Incoming polar angle";
const std::string SphericalCoordinateSystem::ANGLE1_NAME = "Incoming azimuthal angle";
const std::string SphericalCoordinateSystem::ANGLE2_NAME = "Outgoing polar angle";
const std::string SphericalCoordinateSystem::ANGLE3_NAME = "Outgoing azimuthal angle";

const float SphericalCoordinateSystem::MAX_ANGLE0 = PI_2_F      - EPSILON_F * PI_2_F;
const float SphericalCoordinateSystem::MAX_ANGLE1 = 2.0f * PI_F - EPSILON_F * 2.0f * PI_F;
const float SphericalCoordinateSystem::MAX_ANGLE2 = PI_2_F      - EPSILON_F * PI_2_F;
const float SphericalCoordinateSystem::MAX_ANGLE3 = 2.0f * PI_F - EPSILON_F * 2.0f * PI_F;
