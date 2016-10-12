// =================================================================== //
// Copyright (C) 2016 Kimura Ryo                                       //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#include <libbsdf/ReflectanceModel/ReflectanceModel.h>

using namespace lb;

ReflectanceModel::Parameter::Parameter(const std::string& name, float* value)
{
    name_ = name;
    type_ = FLOAT_PARAMETER;
    value_.scalar = value;
}

ReflectanceModel::Parameter::Parameter(const std::string& name, Vec3* value)
{
    name_ = name;
    type_ = VEC3_PARAMETER;
    value_.vec3 = value;
}

Vec3 ReflectanceModel::getBrdfValue(const Vec3& inDir, const Vec3& outDir) const
{
    return getValue(inDir, outDir);
}

std::string ReflectanceModel::getName() const
{
    return "";
}

std::string ReflectanceModel::getDescription() const
{
    return "";
}

ReflectanceModel::~ReflectanceModel()
{
}
