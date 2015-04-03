// =================================================================== //
// Copyright (C) 2015 Kimura Ryo                                       //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#ifndef LIBBSDF_TWO_SIDED_MATERIAL_H
#define LIBBSDF_TWO_SIDED_MATERIAL_H

#include <libbsdf/Brdf/Material.h>

namespace lb {

/*!
 * \class   TwoSidedMaterial
 * \brief   The TwoSidedMaterial class provides the material data of the front and back sides of the surface.
 */
class TwoSidedMaterial
{
public:
    /*!
     * Constructs a two-sided material.
     *
     * \warning \a frontMaterial and \a backMaterial are deleted in destructor.
     */
    TwoSidedMaterial(Material* frontMaterial, Material* backMaterial);

    virtual ~TwoSidedMaterial();

    Material* getFrontMaterial(); /*!< Gets a front side material. */
    Material* getBackMaterial();  /*!< Gets a back side material. */

    const Material* getFrontMaterial() const; /*!< Gets a front side material. */
    const Material* getBackMaterial()  const; /*!< Gets a back side material. */

protected:
    Material* frontMaterial_; /*!< This attribute holds a front side material. */
    Material* backMaterial_;  /*!< This attribute holds a back side material. */

private:
    /*! Copy operator is disabled. */
    TwoSidedMaterial& operator=(const TwoSidedMaterial&);
};

inline Material* TwoSidedMaterial::getFrontMaterial() { return frontMaterial_; }
inline Material* TwoSidedMaterial::getBackMaterial()  { return backMaterial_; }

inline const Material* TwoSidedMaterial::getFrontMaterial() const { return frontMaterial_; }
inline const Material* TwoSidedMaterial::getBackMaterial()  const { return backMaterial_; }

} // namespace lb

#endif // LIBBSDF_TWO_SIDED_MATERIAL_H
