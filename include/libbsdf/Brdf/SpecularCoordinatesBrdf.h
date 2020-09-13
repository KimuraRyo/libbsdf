// =================================================================== //
// Copyright (C) 2014-2020 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#ifndef LIBBSDF_SPECULAR_COORDINATES_BRDF_H
#define LIBBSDF_SPECULAR_COORDINATES_BRDF_H

#include <libbsdf/Brdf/CoordinatesBrdf.h>
#include <libbsdf/Brdf/SphericalCoordinatesBrdf.h>
#include <libbsdf/Common/SpecularCoordinateSystem.h>

namespace lb {

/*!
 * \class   SpecularCoordinatesBrdf
 * \brief   The SpecularCoordinatesBrdf class provides the BRDF of a specular coordinate system.
 *
 * Functions depending on the coordinate system are implemented.
 * \a spec is an abbreviation for specular.
 */
class SpecularCoordinatesBrdf : public CoordinatesBrdf<SpecularCoordinateSystem>
{
private:
    using CoordSys = SpecularCoordinateSystem;
    using BaseBrdf = CoordinatesBrdf<CoordSys>;

public:
    /*!
     * Constructs a spectral BRDF.
     *
     * \param equalIntervalAngles If this parameter is true, angles of sample points are equally-spaced intervals.
     */
    SpecularCoordinatesBrdf(int         numInTheta,
                            int         numInPhi,
                            int         numSpecTheta,
                            int         numSpecPhi,
                            ColorModel  colorModel = RGB_MODEL,
                            int         numWavelengths = 3,
                            bool        equalIntervalAngles = false);

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

    /*!
     * Constructs a BRDF with narrow intervals near specular directions.
     *
     * \param refractiveIndex   Refractive index used for the specular offset of BTDF.
     */
    SpecularCoordinatesBrdf(int         numInTheta,
                            int         numInPhi,
                            int         numSpecTheta,
                            int         numSpecPhi,
                            float       specThetaExponent,
                            ColorModel  colorModel = RGB_MODEL,
                            int         numWavelengths = 3,
                            float       refractiveIndex = 1.0f);

    /*!
     * Constructs a BRDF from lb::SphericalCoordinatesBrdf.
     */
    SpecularCoordinatesBrdf(const SphericalCoordinatesBrdf& brdf,
                            int                             numSpecTheta = 181,
                            int                             numSpecPhi = 73);

    /*! Copies and constructs a BRDF. */
    SpecularCoordinatesBrdf(const SpecularCoordinatesBrdf& brdf);

    virtual ~SpecularCoordinatesBrdf();

    /*! Virtual copy constructor. */
    SpecularCoordinatesBrdf* clone() const override;

    using BaseBrdf::getSpectrum;

    /*! Gets the spectrum of the BRDF at incoming and outgoing directions. */
    Spectrum getSpectrum(const Vec3& inDir, const Vec3& outDir) const override;

    /*! Gets the value of the BRDF at incoming and outgoing directions and the index of wavelength. */
    float getValue(const Vec3& inDir, const Vec3& outDir, int wavelengthIndex) const override;

    /*!
     * Computes incoming and outgoing directions of a Cartesian coordinate system
     * using a set of angle indices.
     */
    void getInOutDirection(int      index0,
                           int      index1,
                           int      index2,
                           int      index3,
                           Vec3*    inDir,
                           Vec3*    outDir) const override;

    /*!
     * Computes outgoing directions of a Cartesian coordinate system
     * using a set of angle indices.
     */
    Vec3 getOutDirection(int index0,
                         int index1,
                         int index2,
                         int index3) const;

    /*! Gets the spectrum of the BRDF at a set of angles. */
    Spectrum getSpectrum(float inTheta,
                         float inPhi,
                         float specTheta,
                         float specPhi);

    /*! Gets the spectrum of the BRDF at a set of angles. */
    Spectrum getSpectrum(float inTheta,
                         float inPhi,
                         float specTheta,
                         float specPhi) const;

    /*! Gets the spectrum of the BRDF at a set of angle indices. */
    Spectrum& getSpectrum(int inThetaIndex,
                          int inPhiIndex,
                          int specThetaIndex,
                          int specPhiIndex);

    /*! Gets the spectrum of the BRDF at a set of angle indices. */
    const Spectrum& getSpectrum(int inThetaIndex,
                                int inPhiIndex,
                                int specThetaIndex,
                                int specPhiIndex) const;

    /*! Gets the spectrum of the isotropic BRDF at a set of angle indices. */
    Spectrum& getSpectrum(int inThetaIndex,
                          int specThetaIndex,
                          int specPhiIndex);

    /*! Gets the spectrum of the isotropic BRDF at a set of angle indices. */
    const Spectrum& getSpectrum(int inThetaIndex,
                                int specThetaIndex,
                                int specPhiIndex) const;

    /*! Sets the spectrum of the BRDF at a set of angle indices. */
    void setSpectrum(int                inThetaIndex,
                     int                inPhiIndex,
                     int                specThetaIndex,
                     int                specPhiIndex,
                     const Spectrum&    spectrum);

    float getInTheta  (int index) const; /*!< Gets the polar angle of an incoming direction. */
    float getInPhi    (int index) const; /*!< Gets the azimuthal angle of an incoming direction. */
    float getSpecTheta(int index) const; /*!< Gets the polar angle of a specular direction. */
    float getSpecPhi  (int index) const; /*!< Gets the azimuthal angle of a specular direction. */

    void setInTheta  (int index, float angle); /*!< Sets the polar angle of an incoming direction. */
    void setInPhi    (int index, float angle); /*!< Sets the azimuthal angle of an incoming direction. */
    void setSpecTheta(int index, float angle); /*!< Sets the polar angle of a specular direction. */
    void setSpecPhi  (int index, float angle); /*!< Sets the azimuthal angle of a specular direction. */

    int getNumInTheta  () const; /*!< Gets the number of polar angles of an incoming direction. */
    int getNumInPhi    () const; /*!< Gets the number of azimuthal angles of an incoming direction. */
    int getNumSpecTheta() const; /*!< Gets the number of polar angles of a specular direction. */
    int getNumSpecPhi  () const; /*!< Gets the number of azimuthal angles of a specular direction. */

    /*! Sets up specular offsets using a refractive index. */
    void setupSpecularOffsets(float refractiveIndex);

    /*! Gets the specular offset at an index. */
    float getSpecularOffset(int index) const;

    /*! Sets the specular offset at an index. */
    void setSpecularOffset(int index, float angle);

    Arrayf&       getSpecularOffsets();       /*!< Gets the array of specular offsets. */
    const Arrayf& getSpecularOffsets() const; /*!< Gets the array of specular offsets. */

    int getNumSpecularOffsets() const; /*!< Gets the number of specular offsets. */

    /*! Gets the specular offset at an incoming polar angle. */
    float getSpecularOffset(float inTheta) const;

    /*!
     * Converts from four angles to incoming and outgoing directions and
     * assigns them to \a inDir and \a outDir.
     * If specular offsets are set, \a inDir and \a outDir are affected by them.
     */
    void toXyz(float inTheta,
               float inPhi,
               float specTheta,
               float specPhi,
               Vec3* inDir,
               Vec3* outDir) const override;

    /*!
     * Converts from incoming and outgoing directions to four angles and assigns them
     * to \a inTheta, \a inPhi, \a specTheta, and \a specPhi.
     * If specular offsets are set, \a specTheta and \a specPhi are affected by them.
     */
    void fromXyz(const Vec3&    inDir,
                 const Vec3&    outDir,
                 float*         inTheta,
                 float*         inPhi,
                 float*         specTheta,
                 float*         specPhi) const override;

    /*!
     * Converts from incoming and outgoing directions to three angles and assigns them
     * to \a inTheta, \a specTheta, and \a specPhi.
     * If specular offsets are set, \a specTheta and \a specPhi are affected by them.
     */
    void fromXyz(const Vec3&    inDir,
                 const Vec3&    outDir,
                 float*         inTheta,
                 float*         specTheta,
                 float*         specPhi) const override;

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
    SpecularCoordinatesBrdf& operator=(const SpecularCoordinatesBrdf&);

    /*!
     * The array of offsets from specular directions.
     * This is efficient to represent the BTDF with refraction.
     * The size of the array is equal to the number of incoming polar angles.
     * If offsets are not used, the array is empty.
     */
    Arrayf specularOffsets_;
};

inline Spectrum SpecularCoordinatesBrdf::getSpectrum(const Vec3& inDir, const Vec3& outDir) const
{
    if (specularOffsets_.size() == 0) {
        return BaseBrdf::getSpectrum(inDir, outDir);
    }

    float inTheta, inPhi, specTheta, specPhi;
    fromXyz(inDir, outDir, &inTheta, &inPhi, &specTheta, &specPhi);

    if (samples_->isIsotropic()) {
        return LinearInterpolator::getSpectrum(*samples_, inTheta, specTheta, specPhi);
    }
    else {
        return LinearInterpolator::getSpectrum(*samples_, inTheta, inPhi, specTheta, specPhi);
    }
}

inline float SpecularCoordinatesBrdf::getValue(const Vec3& inDir, const Vec3& outDir, int wavelengthIndex) const
{
    if (specularOffsets_.size() == 0) {
        return BaseBrdf::getValue(inDir, outDir, wavelengthIndex);
    }

    float inTheta, inPhi, specTheta, specPhi;
    fromXyz(inDir, outDir, &inTheta, &inPhi, &specTheta, &specPhi);

    if (samples_->isIsotropic()) {
        return LinearInterpolator::getValue(*samples_, inTheta, specTheta, specPhi, wavelengthIndex);
    }
    else {
        return LinearInterpolator::getValue(*samples_, inTheta, inPhi, specTheta, specPhi, wavelengthIndex);
    }
}

inline void SpecularCoordinatesBrdf::getInOutDirection(int      index0,
                                                       int      index1,
                                                       int      index2,
                                                       int      index3,
                                                       Vec3*    inDir,
                                                       Vec3*    outDir) const
{
    if (specularOffsets_.size() == 0) {
        BaseBrdf::getInOutDirection(index0, index1, index2, index3, inDir, outDir);
        return;
    }

    toXyz(getInTheta(index0), getInPhi(index1), getSpecTheta(index2), getSpecPhi(index3),
          inDir, outDir);
}

inline Vec3 SpecularCoordinatesBrdf::getOutDirection(int index0,
                                                     int index1,
                                                     int index2,
                                                     int index3) const
{
    float inTheta = getInTheta(index0);
    float inPhi   = getInPhi(index1);

    float offset = getSpecularOffset(inTheta);
    Vec3 outDir = CoordSys::toOutDirXyz(inTheta + offset,
                                        inPhi,
                                        getSpecTheta(index2),
                                        getSpecPhi(index3));
    return outDir;
}

inline Spectrum SpecularCoordinatesBrdf::getSpectrum(float inTheta,
                                                     float inPhi,
                                                     float specTheta,
                                                     float specPhi)
{
    return LinearInterpolator::getSpectrum(*samples_, inTheta, inPhi, specTheta, specPhi);
}

inline Spectrum SpecularCoordinatesBrdf::getSpectrum(float inTheta,
                                                     float inPhi,
                                                     float specTheta,
                                                     float specPhi) const
{
    return LinearInterpolator::getSpectrum(*samples_, inTheta, inPhi, specTheta, specPhi);
}

inline Spectrum& SpecularCoordinatesBrdf::getSpectrum(int inThetaIndex,
                                                      int inPhiIndex,
                                                      int specThetaIndex,
                                                      int specPhiIndex)
{
    return samples_->getSpectrum(inThetaIndex, inPhiIndex, specThetaIndex, specPhiIndex);
}

inline const Spectrum& SpecularCoordinatesBrdf::getSpectrum(int inThetaIndex,
                                                            int inPhiIndex,
                                                            int specThetaIndex,
                                                            int specPhiIndex) const
{
    return samples_->getSpectrum(inThetaIndex, inPhiIndex, specThetaIndex, specPhiIndex);
}

inline Spectrum& SpecularCoordinatesBrdf::getSpectrum(int inThetaIndex,
                                                      int specThetaIndex,
                                                      int specPhiIndex)
{
    return samples_->getSpectrum(inThetaIndex, specThetaIndex, specPhiIndex);
}

inline const Spectrum&  SpecularCoordinatesBrdf::getSpectrum(int inThetaIndex,
                                                             int specThetaIndex,
                                                             int specPhiIndex) const
{
    return samples_->getSpectrum(inThetaIndex, specThetaIndex, specPhiIndex);
}

inline void SpecularCoordinatesBrdf::setSpectrum(int                inThetaIndex,
                                                 int                inPhiIndex,
                                                 int                specThetaIndex,
                                                 int                specPhiIndex,
                                                 const Spectrum&    spectrum)
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

inline int SpecularCoordinatesBrdf::getNumInTheta  () const { return samples_->getNumAngles0(); }
inline int SpecularCoordinatesBrdf::getNumInPhi    () const { return samples_->getNumAngles1(); }
inline int SpecularCoordinatesBrdf::getNumSpecTheta() const { return samples_->getNumAngles2(); }
inline int SpecularCoordinatesBrdf::getNumSpecPhi  () const { return samples_->getNumAngles3(); }

inline float SpecularCoordinatesBrdf::getSpecularOffset(int index) const
{
    return (specularOffsets_.size() == 0) ? 0.0f : specularOffsets_[index];
}

inline void SpecularCoordinatesBrdf::setSpecularOffset(int index, float angle)
{
    if (specularOffsets_.size() == 0) {
        specularOffsets_.resize(getNumInTheta());
    }
    specularOffsets_[index] = angle;
}

inline       Arrayf& SpecularCoordinatesBrdf::getSpecularOffsets()       { return specularOffsets_; }
inline const Arrayf& SpecularCoordinatesBrdf::getSpecularOffsets() const { return specularOffsets_; }

inline int SpecularCoordinatesBrdf::getNumSpecularOffsets() const
{
    return static_cast<int>(specularOffsets_.size());
}

inline float SpecularCoordinatesBrdf::getSpecularOffset(float inTheta) const
{
    if (specularOffsets_.size() == 0) return 0.0f;

    int lIdx0; // index of the lower bound sample point
    int uIdx0; // index of the upper bound sample point
    float lowerAngle0;
    float upperAngle0;

    array_util::findBounds(samples_->getAngles0(), inTheta, samples_->isEqualIntervalAngles0(),
                           &lIdx0, &uIdx0, &lowerAngle0, &upperAngle0);

    float interval = std::max(upperAngle0 - lowerAngle0, EPSILON_F);
    float weight = (inTheta - lowerAngle0) / interval;
    return (specularOffsets_[uIdx0] - specularOffsets_[lIdx0]) * weight + specularOffsets_[lIdx0];
}

inline void SpecularCoordinatesBrdf::toXyz(float inTheta,
                                           float inPhi,
                                           float specTheta,
                                           float specPhi,
                                           Vec3* inDir,
                                           Vec3* outDir) const
{
    *inDir = SphericalCoordinateSystem::toXyz(inTheta, inPhi);

    float offset = getSpecularOffset(inTheta);
    *outDir = CoordSys::toOutDirXyz(inTheta + offset,
                                    inPhi,
                                    specTheta,
                                    specPhi);
}

inline void SpecularCoordinatesBrdf::fromXyz(const Vec3&    inDir,
                                             const Vec3&    outDir,
                                             float*         inTheta,
                                             float*         inPhi,
                                             float*         specTheta,
                                             float*         specPhi) const
{
    SphericalCoordinateSystem::fromXyz(inDir, inTheta, inPhi);

    float offset = getSpecularOffset(*inTheta);
    float offsetInTheta = clamp(*inTheta + offset, 0.0f, CoordSys::MAX_ANGLE0);
    CoordSys::fromOutDirXyz(outDir, offsetInTheta, *inPhi, specTheta, specPhi);
}

inline void SpecularCoordinatesBrdf::fromXyz(const Vec3&    inDir,
                                             const Vec3&    outDir,
                                             float*         inTheta,
                                             float*         specTheta,
                                             float*         specPhi) const
{
    float inPhi;
    SphericalCoordinateSystem::fromXyz(inDir, inTheta, &inPhi);

    float offset = getSpecularOffset(*inTheta);
    float offsetInTheta = clamp(*inTheta + offset, 0.0f, CoordSys::MAX_ANGLE0);
    CoordSys::fromOutDirXyz(outDir, offsetInTheta, inPhi, specTheta, specPhi);
}

} // namespace lb

#endif // LIBBSDF_SPECULAR_COORDINATES_BRDF_H
