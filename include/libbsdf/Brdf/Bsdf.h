// =================================================================== //
// Copyright (C) 2015-2019 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#ifndef LIBBSDF_BSDF_H
#define LIBBSDF_BSDF_H

#include <memory>

#include <libbsdf/Brdf/Brdf.h>
#include <libbsdf/Brdf/Btdf.h>

namespace lb {

/*!
 * \class   Bsdf
 * \brief   The Bsdf class provides the BSDF data structure including BRDF and BTDF.
 */
class Bsdf
{
public:
    /*!
     * Constructs a BSDF.
     *
     * \warning \a brdf and \a btdf are deleted in destructor.
     */
    Bsdf(std::shared_ptr<Brdf> brdf, std::shared_ptr<Btdf> btdf);

    virtual ~Bsdf();

    std::shared_ptr<Brdf> getBrdf(); /*!< Gets the BRDF data. */
    std::shared_ptr<Btdf> getBtdf(); /*!< Gets the BTDF data. */

    const std::shared_ptr<Brdf> getBrdf() const; /*!< Gets the BRDF data. */
    const std::shared_ptr<Btdf> getBtdf() const; /*!< Gets the BTDF data. */

protected:
    /*! This attribute holds the BRDF data including angles, wavelengths, and spectra. */
    std::shared_ptr<Brdf> brdf_;

    /*! This attribute holds the BTDF data including angles, wavelengths, and spectra. */
    std::shared_ptr<Btdf> btdf_;

private:
    /*! Copy operator is disabled. */
    Bsdf& operator=(const Bsdf&);
};

inline std::shared_ptr<Brdf> Bsdf::getBrdf() { return brdf_; }
inline std::shared_ptr<Btdf> Bsdf::getBtdf() { return btdf_; }

inline const std::shared_ptr<Brdf> Bsdf::getBrdf() const { return brdf_; }
inline const std::shared_ptr<Btdf> Bsdf::getBtdf() const { return btdf_; }

} // namespace lb

#endif // LIBBSDF_BSDF_H
