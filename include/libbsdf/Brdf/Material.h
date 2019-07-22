// =================================================================== //
// Copyright (C) 2015 Kimura Ryo                                       //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#ifndef LIBBSDF_MATERIAL_H
#define LIBBSDF_MATERIAL_H

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
    Material(Bsdf*          bsdf,
             SampleSet2D*   specularReflectances,
             SampleSet2D*   specularTransmittances,
             SampleSet2D*   reflectionTis = 0,
             SampleSet2D*   transmissionTis = 0);

    virtual ~Material();

    /*! Gets the BSDF data. */
    Bsdf* getBsdf();

    /*! Gets the BSDF data. */
    const Bsdf* getBsdf() const;

    /*! Gets the array of specular reflectance. */
    SampleSet2D* getSpecularReflectances();

    /*! Gets the array of specular reflectance. */
    const SampleSet2D* getSpecularReflectances() const;

    /*! Gets the array of specular transmittance. */
    SampleSet2D* getSpecularTransmittances();

    /*! Gets the array of specular transmittance. */
    const SampleSet2D* getSpecularTransmittances() const;

    /*! Gets the TIS of reflection. */
    SampleSet2D* getReflectionTis();

    /*! Gets the TIS of reflection. */
    const SampleSet2D* getReflectionTis() const;

    /*! Gets the TIS of transmission. */
    SampleSet2D* getTransmissionTis();

    /*! Gets the TIS of transmission. */
    const SampleSet2D* getTransmissionTis() const;

protected:
    Bsdf* bsdf_; /*!< This attribute holds the BSDF data including angles, wavelengths, and spectra. */

    SampleSet2D* specularReflectances_;     /*!< The array of specular reflectance. */
    SampleSet2D* specularTransmittances_;   /*!< The array of specular transmittance. */

    SampleSet2D* reflectionTis_;    /*!< Total Integrated Scatter (TIS) of reflection. */
    SampleSet2D* transmissionTis_;  /*!< Total Integrated Scatter (TIS) of transmission. */

private:
    /*! Copy operator is disabled. */
    Material& operator=(const Material&);
};

inline       Bsdf* Material::getBsdf()       { return bsdf_; }
inline const Bsdf* Material::getBsdf() const { return bsdf_; }

inline       SampleSet2D* Material::getSpecularReflectances()       { return specularReflectances_; }
inline const SampleSet2D* Material::getSpecularReflectances() const { return specularReflectances_; }

inline       SampleSet2D* Material::getSpecularTransmittances()       { return specularTransmittances_; }
inline const SampleSet2D* Material::getSpecularTransmittances() const { return specularTransmittances_; }

inline       SampleSet2D* Material::getReflectionTis()       { return reflectionTis_; }
inline const SampleSet2D* Material::getReflectionTis() const { return reflectionTis_; }

inline       SampleSet2D* Material::getTransmissionTis()       { return transmissionTis_; }
inline const SampleSet2D* Material::getTransmissionTis() const { return transmissionTis_; }

} // namespace lb

#endif // LIBBSDF_MATERIAL_H
