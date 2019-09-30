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
 * The data includes BSDF, specular reflectances and transmittances, and TIS.
 */
class Material
{
public:
    /*!
     * Constructs a material.
     *
     * \warning \a bsdf, \a specularReflectances, \a specularTransmittances,
     * \a reflectionTis, and \a transmissionTis are deleted in destructor.
     */
    Material(std::shared_ptr<Bsdf>          bsdf,
             std::shared_ptr<SampleSet2D>   specularReflectances,
             std::shared_ptr<SampleSet2D>   specularTransmittances,
             std::shared_ptr<SampleSet2D>   reflectionTis = nullptr,
             std::shared_ptr<SampleSet2D>   transmissionTis = nullptr);

    virtual ~Material();

    /*! Gets the BSDF data. */
    std::shared_ptr<Bsdf> getBsdf();

    /*! Gets the BSDF data. */
    const std::shared_ptr<Bsdf> getBsdf() const;

    /*! Gets the array of specular reflectance. */
    std::shared_ptr<SampleSet2D> getSpecularReflectances();

    /*! Gets the array of specular reflectance. */
    const std::shared_ptr<SampleSet2D> getSpecularReflectances() const;

    /*! Gets the array of specular transmittance. */
    std::shared_ptr<SampleSet2D> getSpecularTransmittances();

    /*! Gets the array of specular transmittance. */
    const std::shared_ptr<SampleSet2D> getSpecularTransmittances() const;

    /*! Gets the TIS of reflection. */
    std::shared_ptr<SampleSet2D> getReflectionTis();

    /*! Gets the TIS of reflection. */
    const std::shared_ptr<SampleSet2D> getReflectionTis() const;

    /*! Gets the TIS of transmission. */
    std::shared_ptr<SampleSet2D> getTransmissionTis();

    /*! Gets the TIS of transmission. */
    const std::shared_ptr<SampleSet2D> getTransmissionTis() const;

protected:
    std::shared_ptr<Bsdf> bsdf_; /*!< This attribute holds the BSDF data including angles, wavelengths, and spectra. */

    std::shared_ptr<SampleSet2D> specularReflectances_;     /*!< The array of specular reflectance. */
    std::shared_ptr<SampleSet2D> specularTransmittances_;   /*!< The array of specular transmittance. */

    std::shared_ptr<SampleSet2D> reflectionTis_;    /*!< Total Integrated Scatter (TIS) of reflection. */
    std::shared_ptr<SampleSet2D> transmissionTis_;  /*!< Total Integrated Scatter (TIS) of transmission. */

private:
    /*! Copy operator is disabled. */
    Material& operator=(const Material&);
};

inline       std::shared_ptr<Bsdf> Material::getBsdf()       { return bsdf_; }
inline const std::shared_ptr<Bsdf> Material::getBsdf() const { return bsdf_; }

inline       std::shared_ptr<SampleSet2D> Material::getSpecularReflectances()       { return specularReflectances_; }
inline const std::shared_ptr<SampleSet2D> Material::getSpecularReflectances() const { return specularReflectances_; }

inline       std::shared_ptr<SampleSet2D> Material::getSpecularTransmittances()       { return specularTransmittances_; }
inline const std::shared_ptr<SampleSet2D> Material::getSpecularTransmittances() const { return specularTransmittances_; }

inline       std::shared_ptr<SampleSet2D> Material::getReflectionTis()       { return reflectionTis_; }
inline const std::shared_ptr<SampleSet2D> Material::getReflectionTis() const { return reflectionTis_; }

inline       std::shared_ptr<SampleSet2D> Material::getTransmissionTis()       { return transmissionTis_; }
inline const std::shared_ptr<SampleSet2D> Material::getTransmissionTis() const { return transmissionTis_; }

} // namespace lb

#endif // LIBBSDF_MATERIAL_H
