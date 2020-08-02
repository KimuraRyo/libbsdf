// =================================================================== //
// Copyright (C) 2015-2020 Kimura Ryo                                  //
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

    /*! Gets the spectrum of the BTDF at incoming and outgoing directions. */
    Spectrum getSpectrum(const Vec3& inDir, const Vec3& outDir) const;

    std::shared_ptr<Brdf> getBrdf(); /*!< Gets the BRDF data. */
    std::shared_ptr<Btdf> getBtdf(); /*!< Gets the BTDF data. */

    std::shared_ptr<const Brdf> getBrdf() const; /*!< Gets the BRDF data. */
    std::shared_ptr<const Btdf> getBtdf() const; /*!< Gets the BTDF data. */

    void setBrdf(std::shared_ptr<Brdf> brdf); /*!< Sets the BRDF data. */
    void setBtdf(std::shared_ptr<Btdf> btdf); /*!< Sets the BTDF data. */

    bool isEmpty() const; /*!< Returns true if BRDF and BTDF do not exist. */

    /*!
     * Validates spectra, angles, wavelengths, and other attributes.
     * False is returned if the data contains one of the following:
     *   - Infinite or NaN spectrum
     *   - Negative spectrum on a visible hemisphere
     *   - Outside, infinite, or NaN angle
     *   - Negative, infinite, or NaN wavelength
     *
     * \param verbose If this parameter is true, all warnings of spectra are output.
     */
    virtual bool validate(bool verbose = false) const;

protected:
    /*! This attribute holds the BRDF data including angles, wavelengths, and spectra. */
    std::shared_ptr<Brdf> brdf_;

    /*! This attribute holds the BTDF data including angles, wavelengths, and spectra. */
    std::shared_ptr<Btdf> btdf_;

private:
    /*! Copy operator is disabled. */
    Bsdf& operator=(const Bsdf&);
};

inline Spectrum Bsdf::getSpectrum(const Vec3& inDir, const Vec3& outDir) const
{
    using std::abs;

    bool reflected = ((inDir[2] >= Vec3::Scalar(0) && outDir[2] >= Vec3::Scalar(0)) ||
                      (inDir[2] <= Vec3::Scalar(0) && outDir[2] <= Vec3::Scalar(0)));
    if (brdf_ && reflected) {
        return brdf_->getSpectrum(Vec3(inDir[0],  inDir[1],  abs(inDir[2])),
                                  Vec3(outDir[0], outDir[1], abs(outDir[2])));
    }
    else if (btdf_ && !reflected) {
        return btdf_->getSpectrum(inDir, outDir);
    }
    else {
        lbWarn << "[Bsdf::getSpectrum] Failed to get a spectrum.";
        return Spectrum(0);
    }
}

inline std::shared_ptr<Brdf> Bsdf::getBrdf() { return brdf_; }
inline std::shared_ptr<Btdf> Bsdf::getBtdf() { return btdf_; }

inline std::shared_ptr<const Brdf> Bsdf::getBrdf() const { return brdf_; }
inline std::shared_ptr<const Btdf> Bsdf::getBtdf() const { return btdf_; }

inline void Bsdf::setBrdf(std::shared_ptr<Brdf> brdf) { brdf_ = brdf; }
inline void Bsdf::setBtdf(std::shared_ptr<Btdf> btdf) { btdf_ = btdf; }

inline bool Bsdf::isEmpty() const { return !(brdf_ || btdf_); }

} // namespace lb

#endif // LIBBSDF_BSDF_H
