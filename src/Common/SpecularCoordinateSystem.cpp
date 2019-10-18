// =================================================================== //
// Copyright (C) 2014-2019 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#include <libbsdf/Common/SpecularCoordinateSystem.h>

using namespace lb;

const char SpecularCoordinateSystem::ANGLE0_NAME[] = "Incoming polar angle";
const char SpecularCoordinateSystem::ANGLE1_NAME[] = "Incoming azimuthal angle";
const char SpecularCoordinateSystem::ANGLE2_NAME[] = "Specular polar angle";
const char SpecularCoordinateSystem::ANGLE3_NAME[] = "Specular azimuthal angle";

const float SpecularCoordinateSystem::MIN_ANGLE0 = 0.0f;
const float SpecularCoordinateSystem::MIN_ANGLE1 = 0.0f;
const float SpecularCoordinateSystem::MIN_ANGLE2 = 0.0f;
const float SpecularCoordinateSystem::MIN_ANGLE3 = 0.0f;

const float SpecularCoordinateSystem::MAX_ANGLE0 = decrease(PI_2_F);
const float SpecularCoordinateSystem::MAX_ANGLE1 = decrease(TAU_F);
const float SpecularCoordinateSystem::MAX_ANGLE2 = decrease(PI_F);
const float SpecularCoordinateSystem::MAX_ANGLE3 = decrease(TAU_F);
