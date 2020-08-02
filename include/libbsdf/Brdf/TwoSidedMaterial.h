// =================================================================== //
// Copyright (C) 2015-2019 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#ifndef LIBBSDF_TWO_SIDED_MATERIAL_H
#define LIBBSDF_TWO_SIDED_MATERIAL_H

#include <memory>

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
    TwoSidedMaterial(std::shared_ptr<Material> frontMaterial,
                     std::shared_ptr<Material> backMaterial);

    virtual ~TwoSidedMaterial();

    std::shared_ptr<Material> getFrontMaterial(); /*!< Gets a front side material. */
    std::shared_ptr<Material> getBackMaterial();  /*!< Gets a back side material. */

    const std::shared_ptr<Material> getFrontMaterial() const; /*!< Gets a front side material. */
    const std::shared_ptr<Material> getBackMaterial()  const; /*!< Gets a back side material. */

protected:
    std::shared_ptr<Material> frontMaterial_; /*!< This attribute holds a front side material. */
    std::shared_ptr<Material> backMaterial_;  /*!< This attribute holds a back side material. */

private:
    /*! Copy operator is disabled. */
    TwoSidedMaterial& operator=(const TwoSidedMaterial&);
};

inline std::shared_ptr<Material> TwoSidedMaterial::getFrontMaterial() { return frontMaterial_; }
inline std::shared_ptr<Material> TwoSidedMaterial::getBackMaterial()  { return backMaterial_; }

inline std::shared_ptr<const Material> TwoSidedMaterial::getFrontMaterial() const { return frontMaterial_; }
inline std::shared_ptr<const Material> TwoSidedMaterial::getBackMaterial()  const { return backMaterial_; }

inline bool TwoSidedMaterial::isEmpty() const
{
    return !((frontMaterial_ && !frontMaterial_->isEmpty()) ||
             (backMaterial_ && !backMaterial_->isEmpty()));
}

} // namespace lb

#endif // LIBBSDF_TWO_SIDED_MATERIAL_H
