// =================================================================== //
// Copyright (C) 2015-2019 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#ifndef LIBBSDF_MATERIAL_H
#define LIBBSDF_MATERIAL_H

#include <memory>

#include <libbsdf/Brdf/Bsdf.h>
#include <libbsdf/Brdf/SampleSet2D.h>

namespace lb {

/*!
 * \class   Material
 * \brief   The Material class provides the material data structure.
 *
 * The data includes a BSDF, specular reflectances, and specular transmittances.
 */
class Material
{
public:
    /*! Constructs a material. */
    explicit Material(std::shared_ptr<Bsdf>         bsdf,
                      std::shared_ptr<SampleSet2D>  specularReflectances = nullptr,
                      std::shared_ptr<SampleSet2D>  specularTransmittances = nullptr);

    virtual ~Material();

    /*! Gets the BSDF data. */
    std::shared_ptr<Bsdf> getBsdf();

    /*! Gets the BSDF data. */
    std::shared_ptr<const Bsdf> getBsdf() const;

    /*! Sets the BSDF data. */
    void setBsdf(std::shared_ptr<Bsdf> bsdf);

    /*! Gets the array of specular reflectance. */
    std::shared_ptr<SampleSet2D> getSpecularReflectances();

    /*! Gets the array of specular reflectance. */
    std::shared_ptr<const SampleSet2D> getSpecularReflectances() const;

    /*! Sets the array of specular reflectance. */
    void setSpecularReflectances(std::shared_ptr<SampleSet2D> specularReflectances);

    /*! Gets the array of specular transmittance. */
    std::shared_ptr<SampleSet2D> getSpecularTransmittances();

    /*! Gets the array of specular transmittance. */
    std::shared_ptr<const SampleSet2D> getSpecularTransmittances() const;

    /*! Sets the array of specular transmittance. */
    void setSpecularTransmittances(std::shared_ptr<SampleSet2D> specularTransmittances);

    /*! Returns true if BRDF, BTDF, specular reflectance, and specular transmittance do not exist. */
    bool isEmpty() const;

protected:
    std::shared_ptr<Bsdf> bsdf_; /*!< This attribute holds the BSDF data including angles, wavelengths, and spectra. */

    std::shared_ptr<SampleSet2D> specularReflectances_;     /*!< The array of specular reflectance. */
    std::shared_ptr<SampleSet2D> specularTransmittances_;   /*!< The array of specular transmittance. */

private:
    /*! Copy operator is disabled. */
    Material& operator=(const Material&);
};

inline std::shared_ptr<      Bsdf> Material::getBsdf()       { return bsdf_; }
inline std::shared_ptr<const Bsdf> Material::getBsdf() const { return bsdf_; }

inline void Material::setBsdf(std::shared_ptr<Bsdf> bsdf) { bsdf_ = bsdf; }

inline std::shared_ptr<      SampleSet2D> Material::getSpecularReflectances()       { return specularReflectances_; }
inline std::shared_ptr<const SampleSet2D> Material::getSpecularReflectances() const { return specularReflectances_; }

inline void Material::setSpecularReflectances(std::shared_ptr<SampleSet2D> specularReflectances)
{
    specularReflectances_ = specularReflectances;
}

inline std::shared_ptr<      SampleSet2D> Material::getSpecularTransmittances()       { return specularTransmittances_; }
inline std::shared_ptr<const SampleSet2D> Material::getSpecularTransmittances() const { return specularTransmittances_; }

inline void Material::setSpecularTransmittances(std::shared_ptr<SampleSet2D> specularTransmittances)
{
    specularTransmittances_ = specularTransmittances;
}

inline bool Material::isEmpty() const
{
    return !((bsdf_ && !bsdf_->isEmpty()) ||
             specularReflectances_ ||
             specularTransmittances_);
}

} // namespace lb

#endif // LIBBSDF_MATERIAL_H
