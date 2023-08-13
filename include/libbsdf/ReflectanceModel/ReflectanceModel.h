// =================================================================== //
// Copyright (C) 2016-2023 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#ifndef LIBBSDF_REFLECTANCE_MODEL_H
#define LIBBSDF_REFLECTANCE_MODEL_H

#include <vector>

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
    /*!
     * \class   Parameter
     * \brief   The Parameter class provides information of a parameter including a name, type, value, etc.
     */
    class Parameter
    {
    public:
        /*! Types of parameter. */
        enum ParameterType
        {
            REAL_PARAMETER,
            VEC3_PARAMETER,
            INT_PARAMETER
        };

        union ValueUnion
        {
            double* real;
            Vec3*   vec3;
            int*    integer;
        };

        Parameter(const std::string& name,
                  double*            value,
                  const double&      minValue = 0.0,
                  const double&      maxValue = 1.0,
                  const std::string& description = "");

        Parameter(const std::string& name,
                  Vec3*              value,
                  const Vec3&        minValue = Vec3::Zero(),
                  const Vec3&        maxValue = Vec3::Ones(),
                  const std::string& description = "");

        Parameter(const std::string& name,
                  int*               value,
                  int                minValue = 0,
                  int                maxValue = 999,
                  const std::string& description = "");

        ~Parameter();
        Parameter(const Parameter& param);
        Parameter& operator=(const Parameter& param);

        const std::string& getName() const;
        ParameterType      getType() const;

        double* getReal() const;
        Vec3*   getVec3() const;
        int*    getInt() const;

        double* getMinReal() const;
        Vec3*   getMinVec3() const;
        int*    getMinInt() const;

        double* getMaxReal() const;
        Vec3*   getMaxVec3() const;
        int*    getMaxInt() const;

        const std::string& getDescription() const;

    private:
        std::string   name_;
        ParameterType type_;
        ValueUnion    value_;
        ValueUnion    minValue_;
        ValueUnion    maxValue_;
        std::string   description_;
    };

    using Parameters = std::vector<Parameter>;

    virtual ~ReflectanceModel();

    /*! Gets a reflected value with incoming and outgoing directions in tangent space. */
    virtual Vec3 getValue(const Vec3& inDir, const Vec3& outDir) const = 0;

    /*!
     * Gets a BRDF value. If a reflectance model is not an analytical BRDF, it
     * should be converted to a BRDF value.
     */
    virtual Vec3 getBrdfValue(const Vec3& inDir, const Vec3& outDir) const;

    /*! Returns true if this reflectance model is isotropic. */
    virtual bool isIsotropic() const = 0;

    /*! Gets the list of parameters for a reflectance model. */
    Parameters& getParameters();

    /*! Gets the list of parameters for a reflectance model. */
    const Parameters& getParameters() const;

    virtual std::string getName() const;

    virtual std::string getDescription() const;

protected:
    Parameters parameters_;
};

inline const std::string& ReflectanceModel::Parameter::getName() const { return name_; }

inline const std::string& ReflectanceModel::Parameter::getDescription() const
{
    return description_;
}

inline ReflectanceModel::Parameter::ParameterType ReflectanceModel::Parameter::getType() const
{
    return type_;
}

inline double* ReflectanceModel::Parameter::getReal() const { return value_.real; }
inline Vec3* ReflectanceModel::Parameter::getVec3() const { return value_.vec3; }
inline int* ReflectanceModel::Parameter::getInt() const { return value_.integer; }

inline double* ReflectanceModel::Parameter::getMinReal() const { return minValue_.real; }
inline Vec3* ReflectanceModel::Parameter::getMinVec3() const { return minValue_.vec3; }
inline int* ReflectanceModel::Parameter::getMinInt() const { return minValue_.integer; }

inline double* ReflectanceModel::Parameter::getMaxReal() const { return maxValue_.real; }
inline Vec3* ReflectanceModel::Parameter::getMaxVec3() const { return maxValue_.vec3; }
inline int* ReflectanceModel::Parameter::getMaxInt() const { return maxValue_.integer; }

inline ReflectanceModel::Parameters& ReflectanceModel::getParameters() { return parameters_; }

inline const ReflectanceModel::Parameters& ReflectanceModel::getParameters() const
{
    return parameters_;
}

} // namespace lb

#endif // LIBBSDF_REFLECTANCE_MODEL_H
