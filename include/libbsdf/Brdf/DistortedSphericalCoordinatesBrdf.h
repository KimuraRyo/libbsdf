// =================================================================== //
// Copyright (C) 2026 Kimura Ryo                                       //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#ifndef LIBBSDF_DISTORTED_SPHERICAL_COORDINATES_BRDF_H
#define LIBBSDF_DISTORTED_SPHERICAL_COORDINATES_BRDF_H

#include <libbsdf/Brdf/CoordinatesBrdf.h>
#include <libbsdf/Brdf/SphericalCoordinatesBrdf.h>
#include <libbsdf/Common/DistortedSphericalCoordinateSystem.h>

namespace lb {

/*!
 * \class   DistortedSphericalCoordinatesBrdf
 * \brief   The DistortedSphericalCoordinatesBrdf class provides the BRDF of a distorted spherical coordinate system.
 *
 * Functions depending on the coordinate system are implemented.
 * \a dist is an abbreviation for distorted.
 */
class DistortedSphericalCoordinatesBrdf : public CoordinatesBrdf<DistortedSphericalCoordinateSystem>
{
private:
    using CoordSys = DistortedSphericalCoordinateSystem;
    using BaseBrdf = CoordinatesBrdf<CoordSys>;

public:
    /*!
     * Constructs a spectral BRDF.
     *
     * \param equalIntervalAngles If this parameter is true, angles of sample points are equally-spaced intervals.
     */
    DistortedSphericalCoordinatesBrdf(int        numInTheta,
                                      int        numInPhi,
                                      int        numDistTheta,
                                      int        numDistPhi,
                                      ColorModel colorModel = RGB_MODEL,
                                      int        numWavelengths = 3,
                                      bool       equalIntervalAngles = false);

    /*! Constructs a BRDF from lb::Brdf and angle lists. */
    DistortedSphericalCoordinatesBrdf(const Brdf&   brdf,
                                      const Arrayd& inThetaAngles,
                                      const Arrayd& inPhiAngles,
                                      const Arrayd& distThetaAngles,
                                      const Arrayd& distPhiAngles);

    /*! Constructs a BRDF from lb::Brdf and the numbers of angles. Angles are equally-spaced intervals. */
    DistortedSphericalCoordinatesBrdf(const Brdf& brdf,
                                      int         numInTheta,
                                      int         numInPhi,
                                      int         numDistTheta,
                                      int         numDistPhi);

    /*!
     * Constructs a BRDF with narrow intervals near specular directions.
     *
     * \param distAngleExponent Exponent for the distribution of distorted angles. The larger this parameter is, the narrower intervals near specular directions are.
     * \param refractiveIndex   Refractive index used for the specular offset of BTDF.
     */
    DistortedSphericalCoordinatesBrdf(int        numInTheta,
                                      int        numInPhi,
                                      int        numDistTheta,
                                      int        numDistPhi,
                                      double     distAngleExponent,
                                      ColorModel colorModel = RGB_MODEL,
                                      int        numWavelengths = 3,
                                      double     refractiveIndex = 1);

    /*!
     * Constructs a BRDF from lb::SphericalCoordinatesBrdf.
     */
    DistortedSphericalCoordinatesBrdf(const SphericalCoordinatesBrdf& brdf,
                                      int                             numDistTheta = 181,
                                      int                             numDistPhi = 73);

    /*! Copies and constructs a BRDF. */
    DistortedSphericalCoordinatesBrdf(const DistortedSphericalCoordinatesBrdf& brdf);

    virtual ~DistortedSphericalCoordinatesBrdf();

    /*! Virtual copy constructor. */
    DistortedSphericalCoordinatesBrdf* clone() const override;

    using BaseBrdf::getSpectrum;

    /*! Gets the spectrum of the BRDF at incoming and outgoing directions. */
    Spectrum getSpectrum(const Vec3& inDir, const Vec3& outDir) const override;

    /*! Gets the value of the BRDF at incoming and outgoing directions and the index of wavelength. */
    float getValue(const Vec3& inDir, const Vec3& outDir, int wavelengthIndex) const override;

    /*!
     * Computes incoming and outgoing directions of a Cartesian coordinate system
     * using a set of angle indices.
     */
    void getInOutDirection(int   index0,
                           int   index1,
                           int   index2,
                           int   index3,
                           Vec3* inDir,
                           Vec3* outDir) const override;

    /*!
     * Computes outgoing directions of a Cartesian coordinate system
     * using a set of angle indices.
     */
    Vec3 getOutDirection(int index0, int index1, int index2, int index3) const;

    /*! Gets the spectrum of the BRDF at a set of angles. */
    Spectrum getSpectrum(double inTheta, double inPhi, double distTheta, double distPhi);

    /*! Gets the spectrum of the BRDF at a set of angles. */
    Spectrum getSpectrum(double inTheta, double inPhi, double distTheta, double distPhi) const;

    /*! Gets the spectrum of the BRDF at a set of angle indices. */
    Spectrum& getSpectrum(int inThetaIndex, int inPhiIndex, int distThetaIndex, int distPhiIndex);

    /*! Gets the spectrum of the BRDF at a set of angle indices. */
    const Spectrum&
    getSpectrum(int inThetaIndex, int inPhiIndex, int distThetaIndex, int distPhiIndex) const;

    /*! Gets the spectrum of the isotropic BRDF at a set of angle indices. */
    Spectrum& getSpectrum(int inThetaIndex, int distThetaIndex, int distPhiIndex);

    /*! Gets the spectrum of the isotropic BRDF at a set of angle indices. */
    const Spectrum& getSpectrum(int inThetaIndex, int distThetaIndex, int distPhiIndex) const;

    /*! Sets the spectrum of the BRDF at a set of angle indices. */
    void setSpectrum(int             inThetaIndex,
                     int             inPhiIndex,
                     int             distThetaIndex,
                     int             distPhiIndex,
                     const Spectrum& spectrum);

    /*! Gets the polar angle of an incoming direction. */
    double getInTheta(int index) const;

    /*! Gets the azimuthal angle of an incoming direction. */
    double getInPhi(int index) const;

    /*! Gets the polar angle of an outgoing direction in a distorted spherical coordinate system. */
    double getDistTheta(int index) const;

    /*! Gets the azimuthal angle of an outgoing direction in a distorted spherical coordinate system. */
    double getDistPhi(int index) const;

    /*! Sets the polar angle of an incoming direction. */
    void setInTheta(int index, double angle);

    /*! Sets the azimuthal angle of an incoming direction. */
    void setInPhi(int index, double angle);

    /*! Sets the polar angle of an outgoing direction in a distorted spherical coordinate system. */
    void setDistTheta(int index, double angle);

    /*! Sets the azimuthal angle of an outgoing direction in a distorted spherical coordinate system. */
    void setDistPhi(int index, double angle);

    /*! Gets the number of polar angles of an incoming direction. */
    int getNumInTheta() const;

    /*! Gets the number of azimuthal angles of an incoming direction. */
    int getNumInPhi() const;

    /*! Gets the number of polar angles of an outgoing direction in a distorted spherical coordinate system. */
    int getNumDistTheta() const;

    /*! Gets the number of azimuthal angles of an outgoing direction in a distorted spherical coordinate system. */
    int getNumDistPhi() const;

    /*! Sets up specular offsets using a refractive index. */
    void setupSpecularOffsets(double refractiveIndex);

    /*! Gets the specular offset at an index. */
    double getSpecularOffset(int index) const;

    /*! Sets the specular offset at an index. */
    void setSpecularOffset(int index, double angle);

    Arrayd&       getSpecularOffsets();       /*!< Gets the array of specular offsets. */
    const Arrayd& getSpecularOffsets() const; /*!< Gets the array of specular offsets. */

    int getNumSpecularOffsets() const; /*!< Gets the number of specular offsets. */

    /*! Gets the specular offset at an incoming polar angle. */
    double getSpecularOffset(double inTheta) const;

    /*!
     * Converts from four angles to incoming and outgoing directions and
     * assigns them to \a inDir and \a outDir.
     * If specular offsets are set, \a inDir and \a outDir are affected by them.
     */
    void toXyz(double inTheta,
               double inPhi,
               double distTheta,
               double distPhi,
               Vec3*  inDir,
               Vec3*  outDir) const override;

    /*!
     * Converts from incoming and outgoing directions to four angles and assigns them
     * to \a inTheta, \a inPhi, \a distTheta, and \a distPhi.
     * If specular offsets are set, \a distTheta and \a distPhi are affected by them.
     */
    void fromXyz(const Vec3& inDir,
                 const Vec3& outDir,
                 double*     inTheta,
                 double*     inPhi,
                 double*     distTheta,
                 double*     distPhi) const override;

    /*!
     * Converts from incoming and outgoing directions to three angles and assigns them
     * to \a inTheta, \a distTheta, and \a distPhi.
     * If specular offsets are set, \a distTheta and \a distPhi are affected by them.
     */
    void fromXyz(const Vec3& inDir,
                 const Vec3& outDir,
                 double*     inTheta,
                 double*     distTheta,
                 double*     distPhi) const override;

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
    bool validate(bool verbose = false) const override;

    /*!
     * Expands minimum angles to MIN_ANGLE and maximum angles to MAX_ANGLE,
     * and constructs the extrapolated sample set.
     */
    bool expandAngles(bool angle0Expanded = true,
                      bool angle1Expanded = true,
                      bool angle2Expanded = true,
                      bool angle3Expanded = true) override;

private:
    /*! Copy operator is disabled. */
    DistortedSphericalCoordinatesBrdf& operator=(const DistortedSphericalCoordinatesBrdf&) = delete;

    /*!
     * The array of offsets from specular directions.
     * This is efficient to represent the BTDF with refraction.
     * The size of the array is equal to the number of incoming polar angles.
     * If offsets are not used, the array is empty.
     */
    Arrayd specularOffsets_;
};

inline Spectrum DistortedSphericalCoordinatesBrdf::getSpectrum(const Vec3& inDir,
                                                               const Vec3& outDir) const
{
    if (specularOffsets_.size() == 0) {
        return BaseBrdf::getSpectrum(inDir, outDir);
    }

    if (samples_->isIsotropic()) {
        double inTheta, distTheta, distPhi;
        fromXyz(inDir, outDir, &inTheta, &distTheta, &distPhi);
        return LinearInterpolator::getSpectrum(*samples_, inTheta, distTheta, distPhi);
    }
    else {
        double inTheta, inPhi, distTheta, distPhi;
        fromXyz(inDir, outDir, &inTheta, &inPhi, &distTheta, &distPhi);
        return LinearInterpolator::getSpectrum(*samples_, inTheta, inPhi, distTheta, distPhi);
    }
}

inline float DistortedSphericalCoordinatesBrdf::getValue(const Vec3& inDir,
                                                         const Vec3& outDir,
                                                         int         wavelengthIndex) const
{
    if (specularOffsets_.size() == 0) {
        return BaseBrdf::getValue(inDir, outDir, wavelengthIndex);
    }

    if (samples_->isIsotropic()) {
        double inTheta, distTheta, distPhi;
        fromXyz(inDir, outDir, &inTheta, &distTheta, &distPhi);
        return LinearInterpolator::getValue(*samples_, inTheta, distTheta, distPhi,
                                            wavelengthIndex);
    }
    else {
        double inTheta, inPhi, distTheta, distPhi;
        fromXyz(inDir, outDir, &inTheta, &inPhi, &distTheta, &distPhi);
        return LinearInterpolator::getValue(*samples_, inTheta, inPhi, distTheta, distPhi,
                                            wavelengthIndex);
    }
}

inline void DistortedSphericalCoordinatesBrdf::getInOutDirection(int   index0,
                                                                 int   index1,
                                                                 int   index2,
                                                                 int   index3,
                                                                 Vec3* inDir,
                                                                 Vec3* outDir) const
{
    if (specularOffsets_.size() == 0) {
        BaseBrdf::getInOutDirection(index0, index1, index2, index3, inDir, outDir);
        return;
    }

    toXyz(getInTheta(index0), getInPhi(index1), getDistTheta(index2), getDistPhi(index3), inDir,
          outDir);
}

inline Vec3 DistortedSphericalCoordinatesBrdf::getOutDirection(int index0,
                                                               int index1,
                                                               int index2,
                                                               int index3) const
{
    double inTheta = getInTheta(index0);
    double offset = getSpecularOffset(inTheta);

    Vec3 inDir, outDir;
    CoordSys::toXyz(inTheta + offset, getInPhi(index1), getDistTheta(index2), getDistPhi(index3),
                    &inDir, &outDir);

    return outDir;
}

inline Spectrum DistortedSphericalCoordinatesBrdf::getSpectrum(double inTheta,
                                                               double inPhi,
                                                               double distTheta,
                                                               double distPhi)
{
    return LinearInterpolator::getSpectrum(*samples_, inTheta, inPhi, distTheta, distPhi);
}

inline Spectrum DistortedSphericalCoordinatesBrdf::getSpectrum(double inTheta,
                                                               double inPhi,
                                                               double distTheta,
                                                               double distPhi) const
{
    return LinearInterpolator::getSpectrum(*samples_, inTheta, inPhi, distTheta, distPhi);
}

inline Spectrum& DistortedSphericalCoordinatesBrdf::getSpectrum(int inThetaIndex,
                                                                int inPhiIndex,
                                                                int distThetaIndex,
                                                                int distPhiIndex)
{
    return samples_->getSpectrum(inThetaIndex, inPhiIndex, distThetaIndex, distPhiIndex);
}

inline const Spectrum& DistortedSphericalCoordinatesBrdf::getSpectrum(int inThetaIndex,
                                                                      int inPhiIndex,
                                                                      int distThetaIndex,
                                                                      int distPhiIndex) const
{
    return samples_->getSpectrum(inThetaIndex, inPhiIndex, distThetaIndex, distPhiIndex);
}

inline Spectrum& DistortedSphericalCoordinatesBrdf::getSpectrum(int inThetaIndex,
                                                                int distThetaIndex,
                                                                int distPhiIndex)
{
    return samples_->getSpectrum(inThetaIndex, distThetaIndex, distPhiIndex);
}

inline const Spectrum& DistortedSphericalCoordinatesBrdf::getSpectrum(int inThetaIndex,
                                                                      int distThetaIndex,
                                                                      int distPhiIndex) const
{
    return samples_->getSpectrum(inThetaIndex, distThetaIndex, distPhiIndex);
}

inline void DistortedSphericalCoordinatesBrdf::setSpectrum(int             inThetaIndex,
                                                           int             inPhiIndex,
                                                           int             distThetaIndex,
                                                           int             distPhiIndex,
                                                           const Spectrum& spectrum)
{
    samples_->setSpectrum(inThetaIndex, inPhiIndex, distThetaIndex, distPhiIndex, spectrum);
}

inline double DistortedSphericalCoordinatesBrdf::getInTheta(int index) const
{
    return samples_->getAngle0(index);
}

inline double DistortedSphericalCoordinatesBrdf::getInPhi(int index) const
{
    return samples_->getAngle1(index);
}

inline double DistortedSphericalCoordinatesBrdf::getDistTheta(int index) const
{
    return samples_->getAngle2(index);
}

inline double DistortedSphericalCoordinatesBrdf::getDistPhi(int index) const
{
    return samples_->getAngle3(index);
}

inline void DistortedSphericalCoordinatesBrdf::setInTheta(int index, double angle)
{
    setAngle0(index, angle);
}

inline void DistortedSphericalCoordinatesBrdf::setInPhi(int index, double angle)
{
    setAngle1(index, angle);
}

inline void DistortedSphericalCoordinatesBrdf::setDistTheta(int index, double angle)
{
    setAngle2(index, angle);
}

inline void DistortedSphericalCoordinatesBrdf::setDistPhi(int index, double angle)
{
    setAngle3(index, angle);
}

inline int DistortedSphericalCoordinatesBrdf::getNumInTheta  () const { return samples_->getNumAngles0(); }
inline int DistortedSphericalCoordinatesBrdf::getNumInPhi    () const { return samples_->getNumAngles1(); }
inline int DistortedSphericalCoordinatesBrdf::getNumDistTheta() const { return samples_->getNumAngles2(); }
inline int DistortedSphericalCoordinatesBrdf::getNumDistPhi  () const { return samples_->getNumAngles3(); }

inline double DistortedSphericalCoordinatesBrdf::getSpecularOffset(int index) const
{
    return (specularOffsets_.size() == 0) ? 0 : specularOffsets_[index];
}

inline void DistortedSphericalCoordinatesBrdf::setSpecularOffset(int index, double angle)
{
    if (specularOffsets_.size() == 0) {
        specularOffsets_.resize(getNumInTheta());
    }
    specularOffsets_[index] = angle;
}

inline       Arrayd& DistortedSphericalCoordinatesBrdf::getSpecularOffsets() { return specularOffsets_; }
inline const Arrayd& DistortedSphericalCoordinatesBrdf::getSpecularOffsets() const
{
    return specularOffsets_;
}

inline int DistortedSphericalCoordinatesBrdf::getNumSpecularOffsets() const
{
    return static_cast<int>(specularOffsets_.size());
}

inline double DistortedSphericalCoordinatesBrdf::getSpecularOffset(double inTheta) const
{
    if (specularOffsets_.size() == 0)
        return 0;

    int    lIdx0; // index of the lower bound sample point
    int    uIdx0; // index of the upper bound sample point
    double lowerAngle0;
    double upperAngle0;

    array_util::findBounds(samples_->getAngles0(), inTheta, samples_->isEqualIntervalAngles0(),
                           &lIdx0, &uIdx0, &lowerAngle0, &upperAngle0);

    double interval = std::max(upperAngle0 - lowerAngle0, EPSILON_D);
    double weight = (inTheta - lowerAngle0) / interval;
    return (specularOffsets_[uIdx0] - specularOffsets_[lIdx0]) * weight + specularOffsets_[lIdx0];
}

inline void DistortedSphericalCoordinatesBrdf::toXyz(double inTheta,
                                                     double inPhi,
                                                     double distTheta,
                                                     double distPhi,
                                                     Vec3*  inDir,
                                                     Vec3*  outDir) const
{
    *inDir = SphericalCoordinateSystem::toXyz(inTheta, inPhi);

    // Compute outDir using the specular offset, but avoid overwriting inDir.
    double offset = getSpecularOffset(inTheta);
    Vec3   tempInDir;
    CoordSys::toXyz(inTheta + offset, inPhi, distTheta, distPhi, &tempInDir, outDir);
}

inline void DistortedSphericalCoordinatesBrdf::fromXyz(const Vec3& inDir,
                                                       const Vec3& outDir,
                                                       double*     inTheta,
                                                       double*     inPhi,
                                                       double*     distTheta,
                                                       double*     distPhi) const
{
    if (specularOffsets_.size() == 0) {
        CoordSys::fromXyz(inDir, outDir, inTheta, inPhi, distTheta, distPhi);
        return;
    }

    SphericalCoordinateSystem::fromXyz(inDir, inTheta, inPhi);

    double offset = getSpecularOffset(*inTheta);
    double offsetInTheta = clamp(*inTheta + offset, 0.0, CoordSys::MAX_ANGLE0);
    Vec3   offsetInDir = SphericalCoordinateSystem::toXyz(offsetInTheta, *inPhi);

    double th, ph;
    CoordSys::fromXyz(offsetInDir, outDir, &th, &ph, distTheta, distPhi);
}

inline void DistortedSphericalCoordinatesBrdf::fromXyz(const Vec3& inDir,
                                                       const Vec3& outDir,
                                                       double*     inTheta,
                                                       double*     distTheta,
                                                       double*     distPhi) const
{
    double inPhi;
    if (specularOffsets_.size() == 0) {
        CoordSys::fromXyz(inDir, outDir, inTheta, &inPhi, distTheta, distPhi);
        return;
    }

    SphericalCoordinateSystem::fromXyz(inDir, inTheta, &inPhi);

    double offset = getSpecularOffset(*inTheta);
    double offsetInTheta = clamp(*inTheta + offset, 0.0, CoordSys::MAX_ANGLE0);
    Vec3   offsetInDir = SphericalCoordinateSystem::toXyz(offsetInTheta, inPhi);

    double th;
    CoordSys::fromXyz(offsetInDir, outDir, &th, distTheta, distPhi);
}

} // namespace lb

#endif // LIBBSDF_DISTORTED_SPHERICAL_COORDINATES_BRDF_H