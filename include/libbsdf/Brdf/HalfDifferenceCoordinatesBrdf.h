// =================================================================== //
// Copyright (C) 2014-2018 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#ifndef LIBBSDF_HALF_DIFFERENCE_COORDINATES_BRDF_H
#define LIBBSDF_HALF_DIFFERENCE_COORDINATES_BRDF_H

#include <libbsdf/Brdf/CoordinatesBrdf.h>
#include <libbsdf/Common/HalfDifferenceCoordinateSystem.h>

namespace lb {

/*!
 * \class   HalfDifferenceCoordinatesBrdf
 * \brief   The HalfDifferenceCoordinatesBrdf class provides the BRDF of a half difference coordinate system.
 *
 * Functions depending on the coordinate system are implemented.
 * \a diff is an abbreviation for difference.
 */
class HalfDifferenceCoordinatesBrdf : public CoordinatesBrdf<HalfDifferenceCoordinateSystem>
{
private:
    typedef HalfDifferenceCoordinateSystem CoordSys;
    typedef CoordinatesBrdf<CoordSys> BaseBrdf;

public:
    /*! Constructs a BRDF. */
    HalfDifferenceCoordinatesBrdf(int           numHalfTheta,
                                  int           numHalfPhi,
                                  int           numDiffTheta,
                                  int           numDiffPhi,
                                  ColorModel    colorModel = RGB_MODEL,
                                  int           numWavelengths = 3,
                                  bool          equalIntervalAngles = false);

    /*! Constructs a BRDF from lb::Brdf and angle lists. */
    HalfDifferenceCoordinatesBrdf(const Brdf&   brdf,
                                  const Arrayf& halfThetaAngles,
                                  const Arrayf& halfPhiAngles,
                                  const Arrayf& diffThetaAngles,
                                  const Arrayf& diffPhiAngles);

    /*! Constructs a BRDF from lb::Brdf and the numbers of angles. Angles are equally-spaced intervals. */
    HalfDifferenceCoordinatesBrdf(const Brdf&   brdf,
                                  int           numHalfTheta,
                                  int           numHalfPhi,
                                  int           numDiffTheta,
                                  int           numDiffPhi);

    /*! Copies and constructs a BRDF. */
    HalfDifferenceCoordinatesBrdf(const HalfDifferenceCoordinatesBrdf& brdf);

    virtual ~HalfDifferenceCoordinatesBrdf();

    /*! Virtual copy constructor. */
    virtual HalfDifferenceCoordinatesBrdf* clone() const;
    
    using BaseBrdf::getSpectrum;

    /*! Gets the spectrum of the BRDF at a set of angles. */
    Spectrum getSpectrum(float halfTheta, float halfPhi,
                         float diffTheta, float diffPhi);

    /*! Gets the spectrum of the BRDF at a set of angles. */
    Spectrum getSpectrum(float halfTheta, float halfPhi,
                         float diffTheta, float diffPhi) const;

    /*! Gets the spectrum of the BRDF at a set of angle indices. */
    Spectrum& getSpectrum(int halfThetaIndex, int halfPhiIndex,
                          int diffThetaIndex, int diffPhiIndex);

    /*! Gets the spectrum of the BRDF at a set of angle indices. */
    const Spectrum& getSpectrum(int halfThetaIndex, int halfPhiIndex,
                                int diffThetaIndex, int diffPhiIndex) const;

    /*! Gets the spectrum of the isotropic BRDF at a set of angle indices. */
    Spectrum& getSpectrum(int halfThetaIndex,
                          int diffThetaIndex, int diffPhiIndex);

    /*! Gets the spectrum of the isotropic BRDF at a set of angle indices. */
    const Spectrum& getSpectrum(int halfThetaIndex,
                                int diffThetaIndex, int diffPhiIndex) const;

    /*! Sets the spectrum of the BRDF at a set of angle indices. */
    void setSpectrum(int halfThetaIndex, int halfPhiIndex,
                     int diffThetaIndex, int diffPhiIndex,
                     const Spectrum& spectrum);

    float getHalfTheta(int index) const; /*!< Gets the polar angle of a halfway vector. */
    float getHalfPhi  (int index) const; /*!< Gets the azimuthal angle of a halfway vector. */
    float getDiffTheta(int index) const; /*!< Gets the polar angle of a difference vector. */
    float getDiffPhi  (int index) const; /*!< Gets the azimuthal angle of a difference vector. */

    void setHalfTheta(int index, float angle); /*!< Sets the polar angle of a halfway vector. */
    void setHalfPhi  (int index, float angle); /*!< Sets the azimuthal angle of a halfway vector. */
    void setDiffTheta(int index, float angle); /*!< Sets the polar angle of a difference vector. */
    void setDiffPhi  (int index, float angle); /*!< Sets the azimuthal angle of a difference vector. */

    int getNumHalfTheta() const; /*!< Gets the number of polar angles of a halfway vector. */
    int getNumHalfPhi()   const; /*!< Gets the number of azimuthal angles of a halfway vector. */
    int getNumDiffTheta() const; /*!< Gets the number of polar angles of a difference vector. */
    int getNumDiffPhi()   const; /*!< Gets the number of azimuthal angles of a difference vector. */

private:
    /*! Copy operator is disabled. */
    HalfDifferenceCoordinatesBrdf& operator=(const HalfDifferenceCoordinatesBrdf&);
};

inline Spectrum HalfDifferenceCoordinatesBrdf::getSpectrum(float halfTheta, float halfPhi,
                                                           float diffTheta, float diffPhi)
{
    Spectrum sp;
    LinearInterpolator::getSpectrum(*samples_, halfTheta, halfPhi, diffTheta, diffPhi, &sp);
    return sp;
}

inline Spectrum HalfDifferenceCoordinatesBrdf::getSpectrum(float halfTheta, float halfPhi,
                                                           float diffTheta, float diffPhi) const
{
    Spectrum sp;
    LinearInterpolator::getSpectrum(*samples_, halfTheta, halfPhi, diffTheta, diffPhi, &sp);
    return sp;
}

inline Spectrum& HalfDifferenceCoordinatesBrdf::getSpectrum(int halfThetaIndex, int halfPhiIndex,
                                                            int diffThetaIndex, int diffPhiIndex)
{
    return samples_->getSpectrum(halfThetaIndex, halfPhiIndex, diffThetaIndex, diffPhiIndex);
}

inline const Spectrum& HalfDifferenceCoordinatesBrdf::getSpectrum(int halfThetaIndex, int halfPhiIndex,
                                                                  int diffThetaIndex, int diffPhiIndex) const
{
    return samples_->getSpectrum(halfThetaIndex, halfPhiIndex, diffThetaIndex, diffPhiIndex);
}

inline Spectrum& HalfDifferenceCoordinatesBrdf::getSpectrum(int halfThetaIndex,
                                                            int diffThetaIndex, int diffPhiIndex)
{
    return samples_->getSpectrum(halfThetaIndex, diffThetaIndex, diffPhiIndex);
}

inline const Spectrum& HalfDifferenceCoordinatesBrdf::getSpectrum(int halfThetaIndex,
                                                                  int diffThetaIndex, int diffPhiIndex) const
{
    return samples_->getSpectrum(halfThetaIndex, diffThetaIndex, diffPhiIndex);
}

inline void HalfDifferenceCoordinatesBrdf::setSpectrum(int halfThetaIndex, int halfPhiIndex,
                                                       int diffThetaIndex, int diffPhiIndex,
                                                       const Spectrum& spectrum)
{
    samples_->setSpectrum(halfThetaIndex, halfPhiIndex, diffThetaIndex, diffPhiIndex, spectrum);
}

inline float HalfDifferenceCoordinatesBrdf::getHalfTheta(int index) const { return samples_->getAngle0(index); }
inline float HalfDifferenceCoordinatesBrdf::getHalfPhi  (int index) const { return samples_->getAngle1(index); }
inline float HalfDifferenceCoordinatesBrdf::getDiffTheta(int index) const { return samples_->getAngle2(index); }
inline float HalfDifferenceCoordinatesBrdf::getDiffPhi  (int index) const { return samples_->getAngle3(index); }

inline void HalfDifferenceCoordinatesBrdf::setHalfTheta(int index, float angle) { setAngle0(index, angle); }
inline void HalfDifferenceCoordinatesBrdf::setHalfPhi  (int index, float angle) { setAngle1(index, angle); }
inline void HalfDifferenceCoordinatesBrdf::setDiffTheta(int index, float angle) { setAngle2(index, angle); }
inline void HalfDifferenceCoordinatesBrdf::setDiffPhi  (int index, float angle) { setAngle3(index, angle); }

inline int HalfDifferenceCoordinatesBrdf::getNumHalfTheta() const { return samples_->getNumAngles0(); }
inline int HalfDifferenceCoordinatesBrdf::getNumHalfPhi()   const { return samples_->getNumAngles1(); }
inline int HalfDifferenceCoordinatesBrdf::getNumDiffTheta() const { return samples_->getNumAngles2(); }
inline int HalfDifferenceCoordinatesBrdf::getNumDiffPhi()   const { return samples_->getNumAngles3(); }

} // namespace lb

#endif // LIBBSDF_HALF_DIFFERENCE_COORDINATES_BRDF_H
