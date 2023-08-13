// =================================================================== //
// Copyright (C) 2014-2023 Kimura Ryo                                  //
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
    using CoordSys = SphericalCoordinateSystem;
    using BaseBrdf = CoordinatesBrdf<CoordSys>;

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
    SphericalCoordinatesBrdf(const Brdf&   brdf,
                             const Arrayd& inThetaAngles,
                             const Arrayd& inPhiAngles,
                             const Arrayd& outThetaAngles,
                             const Arrayd& outPhiAngles);

    /*! Constructs a BRDF from lb::Brdf and the numbers of angles. Angles are equally-spaced intervals. */
    SphericalCoordinatesBrdf(const Brdf& brdf,
                             int         numInTheta,
                             int         numInPhi,
                             int         numOutTheta,
                             int         numOutPhi);

    /*! Copies and constructs a BRDF. */
    SphericalCoordinatesBrdf(const SphericalCoordinatesBrdf& brdf);

    virtual ~SphericalCoordinatesBrdf();

    /*! Virtual copy constructor. */
    SphericalCoordinatesBrdf* clone() const override;

    using BaseBrdf::getSpectrum;

    /*! Gets the spectrum of the BRDF at a set of angles. */
    Spectrum getSpectrum(double inTheta, double inPhi, double outTheta, double outPhi);

    /*! Gets the spectrum of the BRDF at a set of angles. */
    Spectrum getSpectrum(double inTheta, double inPhi, double outTheta, double outPhi) const;

    /*! Gets the spectrum of the BRDF at a set of angle indices. */
    Spectrum& getSpectrum(int inThetaIndex, int inPhiIndex, int outThetaIndex, int outPhiIndex);

    /*! Gets the spectrum of the BRDF at a set of angle indices. */
    const Spectrum&
    getSpectrum(int inThetaIndex, int inPhiIndex, int outThetaIndex, int outPhiIndex) const;

    /*! Gets the spectrum of the isotropic BRDF at a set of angle indices. */
    Spectrum& getSpectrum(int inThetaIndex, int outThetaIndex, int outPhiIndex);

    /*! Gets the spectrum of the isotropic BRDF at a set of angle indices. */
    const Spectrum& getSpectrum(int inThetaIndex, int outThetaIndex, int outPhiIndex) const;

    /*! Sets the spectrum of the BRDF at a set of angle indices. */
    void setSpectrum(int             inThetaIndex,
                     int             inPhiIndex,
                     int             outThetaIndex,
                     int             outPhiIndex,
                     const Spectrum& spectrum);

    double getInTheta(int index) const;  /*!< Gets the polar angle of an incoming direction. */
    double getInPhi(int index) const;    /*!< Gets the azimuthal angle of an incoming direction. */
    double getOutTheta(int index) const; /*!< Gets the polar angle of an outgoing direction. */
    double getOutPhi(int index) const;   /*!< Gets the azimuthal angle of an outgoing direction. */

    void setInTheta(int index, double angle);  /*!< Sets the polar angle of an incoming direction. */
    void setInPhi(int index, double angle);    /*!< Sets the azimuthal angle of an incoming direction. */
    void setOutTheta(int index, double angle); /*!< Sets the polar angle of an outgoing direction. */
    void setOutPhi(int index, double angle);   /*!< Sets the azimuthal angle of an outgoing direction. */

    int getNumInTheta() const;  /*!< Gets the number of polar angles of an incoming direction. */
    int getNumInPhi() const;    /*!< Gets the number of azimuthal angles of an incoming direction. */
    int getNumOutTheta() const; /*!< Gets the number of polar angles of an outgoing direction. */
    int getNumOutPhi() const;   /*!< Gets the number of azimuthal angles of an outgoing direction. */

private:
    /*! Copy operator is disabled. */
    SphericalCoordinatesBrdf& operator=(const SphericalCoordinatesBrdf&);
};

inline Spectrum
SphericalCoordinatesBrdf::getSpectrum(double inTheta, double inPhi, double outTheta, double outPhi)
{
    return LinearInterpolator::getSpectrum(*samples_, inTheta, inPhi, outTheta, outPhi);
}

inline Spectrum SphericalCoordinatesBrdf::getSpectrum(double inTheta,
                                                      double inPhi,
                                                      double outTheta,
                                                      double outPhi) const
{
    return LinearInterpolator::getSpectrum(*samples_, inTheta, inPhi, outTheta, outPhi);
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

inline double SphericalCoordinatesBrdf::getInTheta(int index) const
{
    return samples_->getAngle0(index);
}

inline double SphericalCoordinatesBrdf::getInPhi(int index) const
{
    return samples_->getAngle1(index);
}

inline double SphericalCoordinatesBrdf::getOutTheta(int index) const
{
    return samples_->getAngle2(index);
}

inline double SphericalCoordinatesBrdf::getOutPhi(int index) const
{
    return samples_->getAngle3(index);
}

inline void SphericalCoordinatesBrdf::setInTheta(int index, double angle)
{
    setAngle0(index, angle);
}

inline void SphericalCoordinatesBrdf::setInPhi(int index, double angle)
{
    setAngle1(index, angle);
}

inline void SphericalCoordinatesBrdf::setOutTheta(int index, double angle)
{
    setAngle2(index, angle);
}

inline void SphericalCoordinatesBrdf::setOutPhi(int index, double angle)
{
    setAngle3(index, angle);
}

inline int SphericalCoordinatesBrdf::getNumInTheta () const { return samples_->getNumAngles0(); }
inline int SphericalCoordinatesBrdf::getNumInPhi   () const { return samples_->getNumAngles1(); }
inline int SphericalCoordinatesBrdf::getNumOutTheta() const { return samples_->getNumAngles2(); }
inline int SphericalCoordinatesBrdf::getNumOutPhi  () const { return samples_->getNumAngles3(); }

} // namespace lb

#endif // LIBBSDF_SPHERICAL_COORDINATES_BRDF_H
