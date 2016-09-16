// =================================================================== //
// Copyright (C) 2014-2015 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#include <libbsdf/Common/HalfDifferenceCoordinateSystem.h>

using namespace lb;

const std::string HalfDifferenceCoordinateSystem::ANGLE0_NAME = "Half polar angle";
const std::string HalfDifferenceCoordinateSystem::ANGLE1_NAME = "Half azimuthal angle";
const std::string HalfDifferenceCoordinateSystem::ANGLE2_NAME = "Difference polar angle";
const std::string HalfDifferenceCoordinateSystem::ANGLE3_NAME = "Difference azimuthal angle";

const float HalfDifferenceCoordinateSystem::MIN_ANGLE0 = 0.0f;
const float HalfDifferenceCoordinateSystem::MIN_ANGLE1 = 0.0f;
const float HalfDifferenceCoordinateSystem::MIN_ANGLE2 = 0.0f;
const float HalfDifferenceCoordinateSystem::MIN_ANGLE3 = 0.0f;

const float HalfDifferenceCoordinateSystem::MAX_ANGLE0 = PI_2_F      - EPSILON_F * PI_2_F;
const float HalfDifferenceCoordinateSystem::MAX_ANGLE1 = 2.0f * PI_F - EPSILON_F * 2.0f * PI_F;
const float HalfDifferenceCoordinateSystem::MAX_ANGLE2 = PI_2_F      - EPSILON_F * PI_2_F;
const float HalfDifferenceCoordinateSystem::MAX_ANGLE3 = 2.0f * PI_F - EPSILON_F * 2.0f * PI_F;
