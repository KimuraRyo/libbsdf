// =================================================================== //
// Copyright (C) 2014-2015 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#ifndef LIBBSDF_SPECULAR_COORDINATES_BRDF_H
#define LIBBSDF_SPECULAR_COORDINATES_BRDF_H

#include <libbsdf/Brdf/CoordinatesBrdf.h>
#include <libbsdf/Common/SpecularCoordinateSystem.h>

namespace lb {

/*!
 * \class   SpecularCoordinatesBrdf
 * \brief   The SpecularCoordinatesBrdf class provides the BRDF of a specular coordinate system.
 *
 * Functions depending on the coordinate system are implemented.
 * Spec is an abbreviation for specular.
 */
class SpecularCoordinatesBrdf : public CoordinatesBrdf<SpecularCoordinateSystem>
{
private:
    typedef SpecularCoordinateSystem CoordSys;
    typedef CoordinatesBrdf<CoordSys> BaseBrdf;

public:
    /*! Constructs a spectral BRDF. */
    SpecularCoordinatesBrdf(int                 numInTheta,
                            int                 numInPhi,
                            int                 numSpecTheta,
                            int                 numSpecPhi,
                            ColorModel::Type    colorModel = ColorModel::RGB,
                            int                 numWavelengths = 3,
                            bool                equalIntervalAngles = false);

    /*! Constructs a BRDF from lb::Brdf and angle lists. */
    SpecularCoordinatesBrdf(const Brdf&     brdf,
                            const Arrayf&   inThetaAngles,
                            const Arrayf&   inPhiAngles,
                            const Arrayf&   specThetaAngles,
                            const Arrayf&   specPhiAngles);

    /*! Constructs a BRDF from lb::Brdf and the numbers of angles. Angles are equally-spaced intervals. */
    SpecularCoordinatesBrdf(const Brdf& brdf,
                            int         numInTheta,
                            int         numInPhi,
                            int         numSpecTheta,
                            int         numSpecPhi);

    /*! Copies and constructs a BRDF. */
    SpecularCoordinatesBrdf(const SpecularCoordinatesBrdf& brdf);

    virtual ~SpecularCoordinatesBrdf();
    
    using BaseBrdf::getSpectrum;

    /*! Gets the spectrum of a BRDF at a set of angles. */
    Spectrum getSpectrum(float inTheta, float inPhi,
                         float specTheta, float specPhi) const;

    /*! Gets the spectrum of a BRDF at a set of angle indices. */
    Spectrum& getSpectrum(int inThetaIndex, int inPhiIndex,
                          int specThetaIndex, int specPhiIndex);

    /*! Gets the spectrum of a BRDF at a set of angle indices. */
    const Spectrum& getSpectrum(int inThetaIndex, int inPhiIndex,
                                int specThetaIndex, int specPhiIndex) const;

    /*! Gets the spectrum of an isotropic BRDF at a set of angle indices. */
    Spectrum& getSpectrum(int inThetaIndex,
                          int specThetaIndex, int specPhiIndex);

    /*! Gets the spectrum of an isotropic BRDF at a set of angle indices. */
    const Spectrum& getSpectrum(int inThetaIndex,
                                int specThetaIndex, int specPhiIndex) const;

    /*! Sets the spectrum of a BRDF at a set of angle indices. */
    void setSpectrum(int inThetaIndex, int inPhiIndex,
                     int specThetaIndex, int specPhiIndex,
                     const Spectrum& spectrum);

    float getInTheta  (int index) const; /*!< Gets the polar angle of an incoming direction. */
    float getInPhi    (int index) const; /*!< Gets the azimuthal angle of an incoming direction. */
    float getSpecTheta(int index) const; /*!< Gets the polar angle of a specular direction. */
    float getSpecPhi  (int index) const; /*!< Gets the azimuthal angle of a specular direction. */

    void setInTheta  (int index, float angle); /*!< Sets the polar angle of an incoming direction. */
    void setInPhi    (int index, float angle); /*!< Sets the azimuthal angle of an incoming direction. */
    void setSpecTheta(int index, float angle); /*!< Sets the polar angle of a specular direction. */
    void setSpecPhi  (int index, float angle); /*!< Sets the azimuthal angle of a specular direction. */

    int getNumInTheta()   const; /*!< Gets the number of polar angles of an incoming direction. */
    int getNumInPhi()     const; /*!< Gets the number of azimuthal angles of an incoming direction. */
    int getNumSpecTheta() const; /*!< Gets the number of polar angles of a specular direction. */
    int getNumSpecPhi()   const; /*!< Gets the number of azimuthal angles of a specular direction. */

    /*! Returns true if a BRDF is isotropic. */
    bool isIsotropic() const;

    /*! Fixes the energy conservation of a BRDF with each incoming direction. */
    void fixEnergyConservation();

private:
    /*! Copy operator is disabled. */
    SpecularCoordinatesBrdf& operator=(const SpecularCoordinatesBrdf&);

    /*! Gets the index of spectra from angle indices of an isotropic BRDF. */
    int getIndex(int inThetaIndex,
                 int specThetaIndex, int specPhiIndex) const;
};

inline Spectrum SpecularCoordinatesBrdf::getSpectrum(float inTheta, float inPhi,
                                                     float specTheta, float specPhi) const
{
    Spectrum sp;
    LinearInterpolator::getSpectrum(*samples_, inTheta, inPhi, specTheta, specPhi, &sp);
    return sp;
}

inline Spectrum& SpecularCoordinatesBrdf::getSpectrum(int inThetaIndex, int inPhiIndex,
                                                      int specThetaIndex, int specPhiIndex)
{
    return samples_->getSpectrum(inThetaIndex, inPhiIndex, specThetaIndex, specPhiIndex);
}

inline const Spectrum& SpecularCoordinatesBrdf::getSpectrum(int inThetaIndex, int inPhiIndex,
                                                            int specThetaIndex, int specPhiIndex) const
{
    return samples_->getSpectrum(inThetaIndex, inPhiIndex, specThetaIndex, specPhiIndex);
}

inline Spectrum& SpecularCoordinatesBrdf::getSpectrum(int inThetaIndex,
                                                      int specThetaIndex, int specPhiIndex)
{
    return samples_->getSpectrum(getIndex(inThetaIndex, specThetaIndex, specPhiIndex));
}

inline const Spectrum&  SpecularCoordinatesBrdf::getSpectrum(int inThetaIndex,
                                                             int specThetaIndex, int specPhiIndex) const
{
    return samples_->getSpectrum(getIndex(inThetaIndex, specThetaIndex, specPhiIndex));
}

inline void SpecularCoordinatesBrdf::setSpectrum(int inThetaIndex, int inPhiIndex,
                                                 int specThetaIndex, int specPhiIndex,
                                                 const Spectrum& spectrum)
{
    samples_->setSpectrum(inThetaIndex, inPhiIndex, specThetaIndex, specPhiIndex, spectrum);
}

inline float SpecularCoordinatesBrdf::getInTheta  (int index) const { return samples_->getAngle0(index); }
inline float SpecularCoordinatesBrdf::getInPhi    (int index) const { return samples_->getAngle1(index); }
inline float SpecularCoordinatesBrdf::getSpecTheta(int index) const { return samples_->getAngle2(index); }
inline float SpecularCoordinatesBrdf::getSpecPhi  (int index) const { return samples_->getAngle3(index); }

inline void SpecularCoordinatesBrdf::setInTheta  (int index, float angle) { setAngle0(index, angle); }
inline void SpecularCoordinatesBrdf::setInPhi    (int index, float angle) { setAngle1(index, angle); }
inline void SpecularCoordinatesBrdf::setSpecTheta(int index, float angle) { setAngle2(index, angle); }
inline void SpecularCoordinatesBrdf::setSpecPhi  (int index, float angle) { setAngle3(index, angle); }

inline int SpecularCoordinatesBrdf::getNumInTheta()   const { return samples_->getNumAngles0(); }
inline int SpecularCoordinatesBrdf::getNumInPhi()     const { return samples_->getNumAngles1(); }
inline int SpecularCoordinatesBrdf::getNumSpecTheta() const { return samples_->getNumAngles2(); }
inline int SpecularCoordinatesBrdf::getNumSpecPhi()   const { return samples_->getNumAngles3(); }

inline bool SpecularCoordinatesBrdf::isIsotropic() const { return (getNumInPhi() == 1); }

inline int SpecularCoordinatesBrdf::getIndex(int inThetaIndex,
                                             int specThetaIndex, int specPhiIndex) const
{
    int index = inThetaIndex
              + samples_->getNumAngles0() * specThetaIndex
              + samples_->getNumAngles0() * samples_->getNumAngles2() * specPhiIndex;
    return index;
}

} // namespace lb

#endif // LIBBSDF_SPECULAR_COORDINATES_BRDF_H
