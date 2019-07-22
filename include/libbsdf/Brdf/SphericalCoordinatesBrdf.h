// =================================================================== //
// Copyright (C) 2014-2019 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#ifndef LIBBSDF_SPHERICAL_COORDINATES_BRDF_H
#define LIBBSDF_SPHERICAL_COORDINATES_BRDF_H

#include <libbsdf/Brdf/CoordinatesBrdf.h>
#include <libbsdf/Common/SphericalCoordinateSystem.h>

namespace lb {

/*!
 * \class   SphericalCoordinatesBrdf
 * \brief   The SphericalCoordinatesBrdf class provides the BRDF of a spherical coordinate system.
 *
 * Functions depending on the coordinate system are implemented.
 */
class SphericalCoordinatesBrdf : public CoordinatesBrdf<SphericalCoordinateSystem>
{
private:
    typedef SphericalCoordinateSystem CoordSys;
    typedef CoordinatesBrdf<CoordSys> BaseBrdf;

public:
    /*! Constructs a BRDF. */
    SphericalCoordinatesBrdf(int        numInTheta,
                             int        numInPhi,
                             int        numOutTheta,
                             int        numOutPhi,
                             ColorModel colorModel = RGB_MODEL,
                             int        numWavelengths = 3,
                             bool       equalIntervalAngles = false);

    /*! Constructs a BRDF from lb::Brdf and angle lists. */
    SphericalCoordinatesBrdf(const Brdf&    brdf,
                             const Arrayf&  inThetaAngles,
                             const Arrayf&  inPhiAngles,
                             const Arrayf&  outThetaAngles,
                             const Arrayf&  outPhiAngles);

    /*! Constructs a BRDF from lb::Brdf and the numbers of angles. Angles are equally-spaced intervals. */
    SphericalCoordinatesBrdf(const Brdf&    brdf,
                             int            numInTheta,
                             int            numInPhi,
                             int            numOutTheta,
                             int            numOutPhi);

    /*! Copies and constructs a BRDF. */
    SphericalCoordinatesBrdf(const SphericalCoordinatesBrdf& brdf);

    virtual ~SphericalCoordinatesBrdf();

    /*! Virtual copy constructor. */
    virtual SphericalCoordinatesBrdf* clone() const;
    
    using BaseBrdf::getSpectrum;

    /*! Gets the spectrum of the BRDF at a set of angles. */
    Spectrum getSpectrum(float inTheta,
                         float inPhi,
                         float outTheta,
                         float outPhi);

    /*! Gets the spectrum of the BRDF at a set of angles. */
    Spectrum getSpectrum(float inTheta,
                         float inPhi,
                         float outTheta,
                         float outPhi) const;
    
    /*! Gets the spectrum of the BRDF at a set of angle indices. */
    Spectrum& getSpectrum(int inThetaIndex,
                          int inPhiIndex,
                          int outThetaIndex,
                          int outPhiIndex);

    /*! Gets the spectrum of the BRDF at a set of angle indices. */
    const Spectrum& getSpectrum(int inThetaIndex,
                                int inPhiIndex,
                                int outThetaIndex,
                                int outPhiIndex) const;

    /*! Gets the spectrum of the isotropic BRDF at a set of angle indices. */
    Spectrum& getSpectrum(int inThetaIndex,
                          int outThetaIndex,
                          int outPhiIndex);

    /*! Gets the spectrum of the isotropic BRDF at a set of angle indices. */
    const Spectrum& getSpectrum(int inThetaIndex,
                                int outThetaIndex,
                                int outPhiIndex) const;

    /*! Sets the spectrum of the BRDF at a set of angle indices. */
    void setSpectrum(int                inThetaIndex,
                     int                inPhiIndex,
                     int                outThetaIndex,
                     int                outPhiIndex,
                     const Spectrum&    spectrum);

    float getInTheta (int index) const; /*!< Gets the polar angle of an incoming direction. */
    float getInPhi   (int index) const; /*!< Gets the azimuthal angle of an incoming direction. */
    float getOutTheta(int index) const; /*!< Gets the polar angle of an outgoing direction. */
    float getOutPhi  (int index) const; /*!< Gets the azimuthal angle of an outgoing direction. */

    void setInTheta (int index, float angle); /*!< Sets the polar angle of an incoming direction. */
    void setInPhi   (int index, float angle); /*!< Sets the azimuthal angle of an incoming direction. */
    void setOutTheta(int index, float angle); /*!< Sets the polar angle of an outgoing direction. */
    void setOutPhi  (int index, float angle); /*!< Sets the azimuthal angle of an outgoing direction. */

    int getNumInTheta () const; /*!< Gets the number of polar angles of an incoming direction. */
    int getNumInPhi   () const; /*!< Gets the number of azimuthal angles of an incoming direction. */
    int getNumOutTheta() const; /*!< Gets the number of polar angles of an outgoing direction. */
    int getNumOutPhi  () const; /*!< Gets the number of azimuthal angles of an outgoing direction. */

private:
    /*! Copy operator is disabled. */
    SphericalCoordinatesBrdf& operator=(const SphericalCoordinatesBrdf&);
};

inline Spectrum SphericalCoordinatesBrdf::getSpectrum(float inTheta,
                                                      float inPhi,
                                                      float outTheta,
                                                      float outPhi)
{
    Spectrum sp;
    LinearInterpolator::getSpectrum(*samples_, inTheta, inPhi, outTheta, outPhi, &sp);
    return sp;
}

inline Spectrum SphericalCoordinatesBrdf::getSpectrum(float inTheta,
                                                      float inPhi,
                                                      float outTheta,
                                                      float outPhi) const
{
    Spectrum sp;
    LinearInterpolator::getSpectrum(*samples_, inTheta, inPhi, outTheta, outPhi, &sp);
    return sp;
}

inline Spectrum& SphericalCoordinatesBrdf::getSpectrum(int inThetaIndex,
                                                       int inPhiIndex,
                                                       int outThetaIndex,
                                                       int outPhiIndex)
{
    return samples_->getSpectrum(inThetaIndex, inPhiIndex, outThetaIndex, outPhiIndex);
}

inline const Spectrum& SphericalCoordinatesBrdf::getSpectrum(int inThetaIndex,
                                                             int inPhiIndex,
                                                             int outThetaIndex,
                                                             int outPhiIndex) const
{
    return samples_->getSpectrum(inThetaIndex, inPhiIndex, outThetaIndex, outPhiIndex);
}

inline Spectrum& SphericalCoordinatesBrdf::getSpectrum(int inThetaIndex,
                                                       int outThetaIndex,
                                                       int outPhiIndex)
{
    return samples_->getSpectrum(inThetaIndex, outThetaIndex, outPhiIndex);
}

inline const Spectrum& SphericalCoordinatesBrdf::getSpectrum(int inThetaIndex,
                                                             int outThetaIndex,
                                                             int outPhiIndex) const
{
    return samples_->getSpectrum(inThetaIndex, outThetaIndex, outPhiIndex);
}

inline void SphericalCoordinatesBrdf::setSpectrum(int               inThetaIndex,
                                                  int               inPhiIndex,
                                                  int               outThetaIndex,
                                                  int               outPhiIndex,
                                                  const Spectrum&   spectrum)
{
    samples_->setSpectrum(inThetaIndex, inPhiIndex, outThetaIndex, outPhiIndex, spectrum);
}

inline float SphericalCoordinatesBrdf::getInTheta (int index) const { return samples_->getAngle0(index); }
inline float SphericalCoordinatesBrdf::getInPhi   (int index) const { return samples_->getAngle1(index); }
inline float SphericalCoordinatesBrdf::getOutTheta(int index) const { return samples_->getAngle2(index); }
inline float SphericalCoordinatesBrdf::getOutPhi  (int index) const { return samples_->getAngle3(index); }

inline void SphericalCoordinatesBrdf::setInTheta (int index, float angle) { setAngle0(index, angle); }
inline void SphericalCoordinatesBrdf::setInPhi   (int index, float angle) { setAngle1(index, angle); }
inline void SphericalCoordinatesBrdf::setOutTheta(int index, float angle) { setAngle2(index, angle); }
inline void SphericalCoordinatesBrdf::setOutPhi  (int index, float angle) { setAngle3(index, angle); }

inline int SphericalCoordinatesBrdf::getNumInTheta () const { return samples_->getNumAngles0(); }
inline int SphericalCoordinatesBrdf::getNumInPhi   () const { return samples_->getNumAngles1(); }
inline int SphericalCoordinatesBrdf::getNumOutTheta() const { return samples_->getNumAngles2(); }
inline int SphericalCoordinatesBrdf::getNumOutPhi  () const { return samples_->getNumAngles3(); }

} // namespace lb

#endif // LIBBSDF_SPHERICAL_COORDINATES_BRDF_H
