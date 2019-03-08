// =================================================================== //
// Copyright (C) 2014-2019 Kimura Ryo                                  //
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
 * \a spec is an abbreviation for specular.
 */
class SpecularCoordinatesBrdf : public CoordinatesBrdf<SpecularCoordinateSystem>
{
private:
    typedef SpecularCoordinateSystem CoordSys;
    typedef CoordinatesBrdf<CoordSys> BaseBrdf;

public:
    /*! Constructs a spectral BRDF. */
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

    /*! Copies and constructs a BRDF. */
    SpecularCoordinatesBrdf(const SpecularCoordinatesBrdf& brdf);

    virtual ~SpecularCoordinatesBrdf();

    /*! Virtual copy constructor. */
    virtual SpecularCoordinatesBrdf* clone() const;
    
    using BaseBrdf::getSpectrum;

    /*! Gets the spectrum of the BRDF at incoming and outgoing directions. */
    virtual Spectrum getSpectrum(const Vec3& inDir, const Vec3& outDir) const;

    /*! Gets the value of the BRDF at incoming and outgoing directions and the index of wavelength. */
    virtual float getValue(const Vec3& inDir, const Vec3& outDir, int wavelengthIndex) const;

    /*!
     * Computes incoming and outgoing directions of a Cartesian coordinate system
     * using a set of angle indices.
     */
    virtual void getInOutDirection(int index0, int index1, int index2, int index3,
                                   Vec3* inDir, Vec3* outDir) const;

    /*! Gets the spectrum of the BRDF at a set of angles. */
    Spectrum getSpectrum(float inTheta, float inPhi,
                         float specTheta, float specPhi);

    /*! Gets the spectrum of the BRDF at a set of angles. */
    Spectrum getSpectrum(float inTheta, float inPhi,
                         float specTheta, float specPhi) const;

    /*! Gets the spectrum of the BRDF at a set of angle indices. */
    Spectrum& getSpectrum(int inThetaIndex, int inPhiIndex,
                          int specThetaIndex, int specPhiIndex);

    /*! Gets the spectrum of the BRDF at a set of angle indices. */
    const Spectrum& getSpectrum(int inThetaIndex, int inPhiIndex,
                                int specThetaIndex, int specPhiIndex) const;

    /*! Gets the spectrum of the isotropic BRDF at a set of angle indices. */
    Spectrum& getSpectrum(int inThetaIndex,
                          int specThetaIndex, int specPhiIndex);

    /*! Gets the spectrum of the isotropic BRDF at a set of angle indices. */
    const Spectrum& getSpectrum(int inThetaIndex,
                                int specThetaIndex, int specPhiIndex) const;

    /*! Sets the spectrum of the BRDF at a set of angle indices. */
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

    float getSpecularOffset(int index) const;       /*!< Gets the specular offset at an index. */
    void setSpecularOffset(int index, float angle); /*!< Sets the specular offset at an index. */

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
    void toXyz(float inTheta, float inPhi,
               float specTheta, float specPhi,
               Vec3* inDir, Vec3* outDir) const;

    /*!
     * Converts from incoming and outgoing directions to four angles and assigns them
     * to \a inTheta, \a inPhi, \a specTheta, and \a specPhi.
     * If specular offsets are set, \a specTheta and \a specPhi are affected by them.
     */
    void fromXyz(const Vec3& inDir, const Vec3& outDir,
                 float* inTheta, float* inPhi,
                 float* specTheta, float* specPhi) const;

    /*!
     * Converts from incoming and outgoing directions to three angles and assigns them
     * to \a inTheta, \a specTheta, and \a specPhi.
     * If specular offsets are set, \a specTheta and \a specPhi are affected by them.
     */
    void fromXyz(const Vec3& inDir, const Vec3& outDir,
                 float* inTheta,
                 float* specTheta, float* specPhi) const;

private:
    /*! Copy operator is disabled. */
    SpecularCoordinatesBrdf& operator=(const SpecularCoordinatesBrdf&);

    /*!
     * The array of offsets from specular polar angles for refraction.
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

    Spectrum sp;
    if (samples_->isIsotropic()) {
        LinearInterpolator::getSpectrum(*samples_, inTheta, specTheta, specPhi, &sp);
    }
    else {
        LinearInterpolator::getSpectrum(*samples_, inTheta, inPhi, specTheta, specPhi, &sp);
    }

    return sp;
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

inline void SpecularCoordinatesBrdf::getInOutDirection(int index0, int index1, int index2, int index3,
                                                       Vec3* inDir, Vec3* outDir) const
{
    if (specularOffsets_.size() == 0) {
        BaseBrdf::getInOutDirection(index0, index1, index2, index3, inDir, outDir);
        return;
    }

    float inTheta = samples_->getAngle0(index0);

    *inDir = SphericalCoordinateSystem::toXyz(inTheta,
                                              samples_->getAngle1(index1));
    inDir->normalize();

    float offset = getSpecularOffset(inTheta);
    *outDir = SpecularCoordinateSystem::toOutDirXyz(inTheta + offset,
                                                    samples_->getAngle1(index1),
                                                    samples_->getAngle2(index2),
                                                    samples_->getAngle3(index3));
    outDir->normalize();
}

inline Spectrum SpecularCoordinatesBrdf::getSpectrum(float inTheta, float inPhi,
                                                     float specTheta, float specPhi)
{
    Spectrum sp;
    LinearInterpolator::getSpectrum(*samples_, inTheta, inPhi, specTheta, specPhi, &sp);
    return sp;
}

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
    return samples_->getSpectrum(inThetaIndex, specThetaIndex, specPhiIndex);
}

inline const Spectrum&  SpecularCoordinatesBrdf::getSpectrum(int inThetaIndex,
                                                             int specThetaIndex, int specPhiIndex) const
{
    return samples_->getSpectrum(inThetaIndex, specThetaIndex, specPhiIndex);
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

inline Arrayf&       SpecularCoordinatesBrdf::getSpecularOffsets()       { return specularOffsets_; }
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

    LinearInterpolator::findBounds(samples_->getAngles0(), inTheta, samples_->isEqualIntervalAngles0(),
                                   &lIdx0, &uIdx0, &lowerAngle0, &upperAngle0);

    float interval = std::max(upperAngle0 - lowerAngle0, EPSILON_F);
    float weight = (inTheta - lowerAngle0) / interval;
    return (specularOffsets_[uIdx0] - specularOffsets_[lIdx0]) * weight + specularOffsets_[lIdx0];
}

inline void SpecularCoordinatesBrdf::toXyz(float inTheta, float inPhi,
                                           float specTheta, float specPhi,
                                           Vec3* inDir, Vec3* outDir) const
{
    float offset = getSpecularOffset(inTheta);

    SpecularCoordinateSystem::toXyz(inTheta + offset, inPhi, specTheta, specPhi, inDir, outDir);
}

inline void SpecularCoordinatesBrdf::fromXyz(const Vec3& inDir, const Vec3& outDir,
                                             float* inTheta, float* inPhi,
                                             float* specTheta, float* specPhi) const
{
    SphericalCoordinateSystem::fromXyz(inDir, inTheta, inPhi);

    float offset = getSpecularOffset(*inTheta);
    SpecularCoordinateSystem::fromOutDirXyz(outDir, *inTheta + offset, *inPhi, specTheta, specPhi);
}

inline void SpecularCoordinatesBrdf::fromXyz(const Vec3& inDir, const Vec3& outDir,
                                             float* inTheta,
                                             float* specTheta, float* specPhi) const
{
    float inPhi;
    SphericalCoordinateSystem::fromXyz(inDir, inTheta, &inPhi);

    float offset = getSpecularOffset(*inTheta);
    SpecularCoordinateSystem::fromOutDirXyz(outDir, *inTheta + offset, inPhi, specTheta, specPhi);
}

} // namespace lb

#endif // LIBBSDF_SPECULAR_COORDINATES_BRDF_H
