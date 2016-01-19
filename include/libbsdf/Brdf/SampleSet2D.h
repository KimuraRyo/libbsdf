// =================================================================== //
// Copyright (C) 2014-2015 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#ifndef LIBBSDF_SAMPLE_SET_2D_H
#define LIBBSDF_SAMPLE_SET_2D_H

#include <libbsdf/Brdf/LinearInterpolator.h>
#include <libbsdf/Brdf/Sampler.h>
#include <libbsdf/Common/Global.h>
#include <libbsdf/Common/Vector.h>

namespace lb {

/*!
 * \class   SampleSet2D
 * \brief   The SampleSet2D provides a 2D sample array.
 *
 * A spherical coordinate system is used.
 * The arrays of angles and wavelengths must be sorted in ascending order.
 */
class SampleSet2D
{
public:
    /*! Constructs a 2D sample array. */
    SampleSet2D(int         numTheta,
                int         numPhi,
                ColorModel  colorModel = RGB_MODEL,
                int         numWavelengths = 3,
                bool        equalIntervalAngles = false);

    /*! Gets the spectrum at a direction. */
    Spectrum getSpectrum(const Vec3& inDir) const;

    /*! Gets the spectrum at a set of angles. */
    Spectrum getSpectrum(float theta, float phi) const;

    /*! Gets the spectrum at a set of angles. */
    Spectrum getSpectrum(float theta, float phi);

    /*! Gets the spectrum at a polar angle. */
    Spectrum getSpectrum(float theta) const;

    /*! Gets the spectrum at a polar angle. */
    Spectrum getSpectrum(float theta);

    /*! Gets the spectrum at a set of angle indices. */
    Spectrum& getSpectrum(int thetaIndex, int phiIndex);

    /*! Gets the spectrum at a set of angle indices. */
    const Spectrum& getSpectrum(int thetaIndex, int phiIndex) const;

    /*! Gets the spectrum at a angle indice. */
    Spectrum& getSpectrum(int thetaIndex);

    /*! Gets the spectrum at a angle indice. */
    const Spectrum& getSpectrum(int thetaIndex) const;

    /*! Sets the spectrum at a set of angle indices. */
    void setSpectrum(int thetaIndex, int phiIndex, const Spectrum& spectrum);

    float getTheta(int index) const; /*!< Gets the polar angle at an index. */
    float getPhi  (int index) const; /*!< Gets the azimuthal angle at an index. */

    void setTheta(int index, float angle); /*!< Sets the polar angle at an index. */
    void setPhi  (int index, float angle); /*!< Sets the azimuthal angle at an index. */

    Arrayf& getThetaArray(); /*!< Gets the array of polar angles. */
    Arrayf& getPhiArray();   /*!< Gets the array of azimuthal angles. */

    const Arrayf& getThetaArray() const; /*!< Gets the array of polar angles. */
    const Arrayf& getPhiArray()   const; /*!< Gets the array of azimuthal angles. */

    int getNumTheta() const; /*!< Gets the number of polar angles. */
    int getNumPhi()   const; /*!< Gets the number of azimuthal angles. */

    /*! Returns true if polar angles are set at equal intervals. */
    bool isEqualIntervalTheta() const;

    /*! Returns true if azimuthal angles are set at equal intervals. */
    bool isEqualIntervalPhi() const;

    /*! Gets the color model. */
    ColorModel getColorModel() const;

    /*! Sets the color model. */
    void setColorModel(ColorModel colorModel);

    /*! Gets the wavelength at an index. */
    float getWavelength(int index) const;

    /*! Sets the wavelength at an index. */
    void setWavelength(int index, float wavelength);

    /*! Gets the array of wavelengths. */
    Arrayf& getWavelengths();

    /*! Gets the array of wavelengths. */
    const Arrayf& getWavelengths() const;

    /*! Gets the number of wavelengths. */
    int getNumWavelengths() const;

    /*! Returns true if the data is isotropic. */
    bool isIsotropic() const;

    /*! Updates angle attributes. */
    void updateAngleAttributes();

    /*! Clamps all angles to minimum and maximum values. */
    void clampAngles();

private:
    SpectrumList spectra_; /*!< The list of spectrum for each direction. */

    Arrayf thetaAngles_; /*!< The array of polar angles. */
    Arrayf phiAngles_;   /*!< The array of azimuthal angles. */

    int numTheta_; /*!< The number of polar angles. */
    int numPhi_;   /*!< The number of azimuthal angles. */

    bool equalIntervalTheta_; /*!< This attribute holds whether polar angles are set at equal intervals. */
    bool equalIntervalPhi_;   /*!< This attribute holds whether azimuthal angles are set at equal intervals. */

    ColorModel colorModel_; /*!< The color model of spectra. */

    Arrayf wavelengths_; /*!< The array of wavelengths. */
};

inline Spectrum SampleSet2D::getSpectrum(const Vec3& inDir) const
{
    Spectrum sp;
    Sampler::getSpectrum<LinearInterpolator>(*this, inDir, &sp);
    return sp;
}

inline Spectrum SampleSet2D::getSpectrum(float theta, float phi) const
{
    Spectrum sp;
    LinearInterpolator::getSpectrum(*this, theta, phi, &sp);
    return sp;
}

inline Spectrum SampleSet2D::getSpectrum(float theta, float phi)
{
    Spectrum sp;
    LinearInterpolator::getSpectrum(*this, theta, phi, &sp);
    return sp;
}

inline Spectrum SampleSet2D::getSpectrum(float theta) const
{
    Spectrum sp;
    LinearInterpolator::getSpectrum(*this, theta, &sp);
    return sp;
}

inline Spectrum SampleSet2D::getSpectrum(float theta)
{
    Spectrum sp;
    LinearInterpolator::getSpectrum(*this, theta, &sp);
    return sp;
}

inline Spectrum& SampleSet2D::getSpectrum(int thetaIndex, int phiIndex)
{
    return spectra_.at(thetaIndex + numTheta_ * phiIndex);
}

inline const Spectrum& SampleSet2D::getSpectrum(int thetaIndex, int phiIndex) const
{
    return spectra_.at(thetaIndex + numTheta_ * phiIndex);
}

inline Spectrum& SampleSet2D::getSpectrum(int thetaIndex)
{
    return spectra_.at(thetaIndex);
}

inline const Spectrum& SampleSet2D::getSpectrum(int thetaIndex) const
{
    return spectra_.at(thetaIndex);
}

inline void SampleSet2D::setSpectrum(int thetaIndex, int phiIndex, const Spectrum& spectrum)
{
    spectra_.at(thetaIndex + numTheta_ * phiIndex) = spectrum;
}

inline float SampleSet2D::getTheta(int index) const { return thetaAngles_[index]; }
inline float SampleSet2D::getPhi(int index)   const { return phiAngles_[index]; }

inline void SampleSet2D::setTheta(int index, float angle)
{
    thetaAngles_[index] = clamp(angle, 0.0f, SphericalCoordinateSystem::MAX_ANGLE0);
    equalIntervalTheta_ = isEqualInterval(thetaAngles_);
}

inline void SampleSet2D::setPhi(int index, float angle)
{
    phiAngles_[index] = clamp(angle, 0.0f, SphericalCoordinateSystem::MAX_ANGLE1);
    equalIntervalPhi_ = isEqualInterval(phiAngles_);
}

inline Arrayf& SampleSet2D::getThetaArray() { return thetaAngles_; }
inline Arrayf& SampleSet2D::getPhiArray()   { return phiAngles_; }

inline const Arrayf& SampleSet2D::getThetaArray() const { return thetaAngles_; }
inline const Arrayf& SampleSet2D::getPhiArray()   const { return phiAngles_; }

inline int SampleSet2D::getNumTheta() const { return numTheta_; }
inline int SampleSet2D::getNumPhi()   const { return numPhi_; }

inline bool SampleSet2D::isEqualIntervalTheta() const { return equalIntervalTheta_; }
inline bool SampleSet2D::isEqualIntervalPhi()   const { return equalIntervalPhi_; }

inline ColorModel SampleSet2D::getColorModel() const { return colorModel_; }

inline void SampleSet2D::setColorModel(ColorModel colorModel)
{
    colorModel_ = colorModel;
}

inline float SampleSet2D::getWavelength(int index) const
{
    return wavelengths_[index];
}

inline void SampleSet2D::setWavelength(int index, float wavelength)
{
    wavelengths_[index] = wavelength;
}

inline Arrayf& SampleSet2D::getWavelengths() { return wavelengths_; }

inline const Arrayf& SampleSet2D::getWavelengths() const { return wavelengths_; }

inline int SampleSet2D::getNumWavelengths() const { return wavelengths_.size(); }

inline bool SampleSet2D::isIsotropic() const { return (numPhi_ == 1); }

} // namespace lb

#endif // LIBBSDF_SAMPLE_SET_2D_H
