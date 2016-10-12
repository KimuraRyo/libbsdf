// =================================================================== //
// Copyright (C) 2016 Kimura Ryo                                       //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#ifndef LIBBSDF_REFLECTANCE_MODEL_H
#define LIBBSDF_REFLECTANCE_MODEL_H

#include<vector>

#include <libbsdf/Common/Vector.h>

namespace lb {

/*!
 * \class   ReflectanceModel
 * \brief   The ReflectanceModel class is the base class of reflectance models.
 *
 * \warning There are reflectance models that are not designed as analytical BRDFs.
 *          The incident cosine law should not be applied for some empirical model (e.g., Phong).
 *          The radiance factor is sometimes used instead of the BRDF in remote sensing (e.g., Minnaert).
 */
class ReflectanceModel
{
public:
    class Parameter
    {
    public:
        enum ParameterType
        {
            FLOAT_PARAMETER,
            VEC3_PARAMETER
        };

        union ValueUnion
        {
            float*  scalar;
            Vec3*   vec3;
        };

        Parameter(const std::string& name, float* value);
        Parameter(const std::string& name, Vec3* value);

        const std::string&  getName() const;
        ParameterType       getType() const;

        float*  getFloat();
        Vec3*   getVec3();

    private:
        std::string     name_;
        ParameterType   type_;
        ValueUnion      value_;
    };

    typedef std::vector<Parameter> Parameters;

    virtual ~ReflectanceModel();
    
    /*! Gets a reflected value with incoming and outgoing directions in tangent space. */
    virtual Vec3 getValue(const Vec3& inDir, const Vec3& outDir) const = 0;

    /*!
     * Gets a BRDF value. If a reflectance model is not an analytical BRDF, it
     * should be converted to a BRDF value.
     */
    virtual Vec3 getBrdfValue(const Vec3& inDir, const Vec3& outDir) const;

    /*! Returns ture if this reflectance model is isotropic. */
    virtual bool isIsotropic() const = 0;

    /*! Gets the list of parameters for a reflectance model. */
    Parameters& getParameters();

    virtual std::string getName() const;

    virtual std::string getDescription() const;

protected:
    Parameters parameters_;
};

inline const std::string& ReflectanceModel::Parameter::getName() const
{
    return name_;
}

inline ReflectanceModel::Parameter::ParameterType ReflectanceModel::Parameter::getType() const
{
    return type_;
}

inline float* ReflectanceModel::Parameter::getFloat()
{
    return value_.scalar;
}

inline Vec3* ReflectanceModel::Parameter::getVec3()
{
    return value_.vec3;
}

inline ReflectanceModel::Parameters& ReflectanceModel::getParameters()
{
    return parameters_;
}

} // namespace lb

#endif // LIBBSDF_REFLECTANCE_MODEL_H
