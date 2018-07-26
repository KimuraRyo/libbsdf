// =================================================================== //
// Copyright (C) 2016-2018 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#include <libbsdf/ReflectanceModel/ReflectanceModel.h>

using namespace lb;

ReflectanceModel::Parameter::Parameter(const std::string&   name,
                                       float*               value,
                                       const float&         minValue,
                                       const float&         maxValue,
                                       const std::string&   description)
{
    name_ = name;
    type_ = FLOAT_PARAMETER;
    value_.real = value;
    minValue_.real = new float(minValue);
    maxValue_.real = new float(maxValue);
    description_ = description;
}

ReflectanceModel::Parameter::Parameter(const std::string&   name,
                                       Vec3*                value,
                                       const Vec3&          minValue,
                                       const Vec3&          maxValue,
                                       const std::string&   description)
{
    name_ = name;
    type_ = VEC3_PARAMETER;
    value_.vec3 = value;
    minValue_.vec3 = new Vec3(minValue);
    maxValue_.vec3 = new Vec3(maxValue);
    description_ = description;
}

ReflectanceModel::Parameter::Parameter(const std::string&   name,
                                       int*                 value,
                                       int                  minValue,
                                       int                  maxValue,
                                       const std::string&   description)
{
    name_ = name;
    type_ = INT_PARAMETER;
    value_.integer = value;
    minValue_.integer = new int(minValue);
    maxValue_.integer = new int(maxValue);
    description_ = description;
}

ReflectanceModel::Parameter::~Parameter()
{
    switch (type_) {
        case FLOAT_PARAMETER:
            delete minValue_.real;
            delete maxValue_.real;
            break;
        case VEC3_PARAMETER:
            delete minValue_.vec3;
            delete maxValue_.vec3;
            break;
        case INT_PARAMETER:
            delete minValue_.integer;
            delete maxValue_.integer;
            break;
        default:
            break;
    }
}

ReflectanceModel::Parameter::Parameter(const ReflectanceModel::Parameter& param)
{
    *this = param;
}

ReflectanceModel::Parameter& ReflectanceModel::Parameter::operator=(const Parameter& param)
{
    if (&param == this) {
        return *this;
    }

    name_ = param.name_;
    type_ = param.type_;
    value_ = param.value_;
    description_ = param.description_;

    switch (type_) {
        case FLOAT_PARAMETER:
            minValue_.real = new float(*param.minValue_.real);
            maxValue_.real = new float(*param.maxValue_.real);
            break;
        case VEC3_PARAMETER:
            minValue_.vec3 = new Vec3(*param.minValue_.vec3);
            maxValue_.vec3 = new Vec3(*param.maxValue_.vec3);
            break;
        case INT_PARAMETER:
            minValue_.integer = new int(*param.minValue_.integer);
            maxValue_.integer = new int(*param.maxValue_.integer);
            break;
        default:
            break;
    }

    return *this;
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
