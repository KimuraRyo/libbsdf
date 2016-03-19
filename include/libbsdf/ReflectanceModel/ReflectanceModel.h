// =================================================================== //
// Copyright (C) 2016 Kimura Ryo                                       //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#ifndef LIBBSDF_REFLECTANCE_MODEL_H
#define LIBBSDF_REFLECTANCE_MODEL_H

#include<map>

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
    typedef std::map<std::string, float*> Parameters;

    virtual ~ReflectanceModel() {}
    
    /*! Gets a reflected value with incoming and outgoing directions in tangent space. */
    virtual float getValue(const Vec3& inDir, const Vec3& outDir) const = 0;

    /*!
     * Gets a BRDF value. If a reflectance model is not an analytical BRDF, it
     * should be converted to a BRDF value.
     */
    virtual float getBrdfValue(const Vec3& inDir, const Vec3& outDir) const
    {
        return getValue(inDir, outDir);
    }

    /*! Returns ture if this reflectance model is isotropic. */
    virtual bool isIsotropic() const = 0;

    /*! Gets the list of parameters for a reflectance model. */
    Parameters& getParameters() { return parameters_; }

    virtual std::string getName() const { return ""; }

    virtual std::string getDescription() const { return ""; }

protected:
    Parameters parameters_;
};

} // namespace lb

#endif // LIBBSDF_REFLECTANCE_MODEL_H
