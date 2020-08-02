// =================================================================== //
// Copyright (C) 2014-2019 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#ifndef LIBBSDF_SAMPLE_SET_2D_H
#define LIBBSDF_SAMPLE_SET_2D_H

#include <libbsdf/Brdf/LinearInterpolator.h>

#include <libbsdf/Common/SphericalCoordinateSystem.h>

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
    /*!
     * Constructs a 2D sample array.
     * The initialization of angles, wavelengths, and spectra is required.
     */
    SampleSet2D(int         numTheta,
                int         numPhi,
                ColorModel  colorModel = RGB_MODEL,
                int         numWavelengths = 3,
                bool        equalIntervalAngles = false);

    ~SampleSet2D();

    /*! Gets the spectrum at a direction. */
    Spectrum getSpectrum(const Vec3& dir) const;

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

    /*! Gets all spectra. */
    SpectrumList& getSpectra();

    /*! Gets all spectra. */
    const SpectrumList& getSpectra() const;

    /*! Gets a direction at a set of angle indices. */
    Vec3 getDirection(int thetaIndex, int phiIndex) const;

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

    /*! Gets the source type. */
    SourceType getSourceType() const;

    /*! Sets the source type. */
    void setSourceType(SourceType type);

    /*! Returns true if the data is isotropic. */
    bool isIsotropic() const;

    /*! \brief Updates angle attributes.
     *
     * Updates the attributes whether angles are set at equal intervals.
     */
    void updateAngleAttributes();

    /*! Resizes the number of angles. Angles and spectra must be initialized. */
    void resizeAngles(int numTheta, int numPhi);

    /*! Resizes the number of wavelengths. Wavelengths and spectra must be initialized. */
    void resizeWavelengths(int numWavelengths);

    /*! Clamps all angles to minimum and maximum values. */
    void clampAngles();

private:
    SpectrumList spectra_; /*!< The list of spectrum for each direction. */

    Arrayf thetaAngles_; /*!< The array of polar angles. */
    Arrayf phiAngles_;   /*!< The array of azimuthal angles. */

    bool equalIntervalTheta_; /*!< This attribute holds whether polar angles are set at equal intervals. */
    bool equalIntervalPhi_;   /*!< This attribute holds whether azimuthal angles are set at equal intervals. */

    ColorModel colorModel_; /*!< The color model of spectra. */

    Arrayf wavelengths_; /*!< The array of wavelengths. */

    SourceType sourceType_; /*!< The data type of source. */
};

inline Spectrum SampleSet2D::getSpectrum(float theta, float phi) const
{
    return LinearInterpolator::getSpectrum(*this, theta, phi);
}

inline Spectrum SampleSet2D::getSpectrum(float theta, float phi)
{
    return LinearInterpolator::getSpectrum(*this, theta, phi);
}

inline Spectrum SampleSet2D::getSpectrum(float theta) const
{
    return LinearInterpolator::getSpectrum(*this, theta);
}

inline Spectrum SampleSet2D::getSpectrum(float theta)
{
    return LinearInterpolator::getSpectrum(*this, theta);
}

inline Spectrum& SampleSet2D::getSpectrum(int thetaIndex, int phiIndex)
{
    return spectra_.at(thetaIndex + thetaAngles_.size() * phiIndex);
}

inline const Spectrum& SampleSet2D::getSpectrum(int thetaIndex, int phiIndex) const
{
    return spectra_.at(thetaIndex + thetaAngles_.size() * phiIndex);
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
    spectra_.at(thetaIndex + thetaAngles_.size() * phiIndex) = spectrum;
}

inline       SpectrumList& SampleSet2D::getSpectra()       { return spectra_; }
inline const SpectrumList& SampleSet2D::getSpectra() const { return spectra_; }

inline Vec3 SampleSet2D::getDirection(int thetaIndex, int phiIndex) const
{
    return SphericalCoordinateSystem::toXyz(thetaAngles_[thetaIndex], phiAngles_[phiIndex]);
}

inline float SampleSet2D::getTheta(int index) const { return thetaAngles_[index]; }
inline float SampleSet2D::getPhi(int index)   const { return phiAngles_[index]; }

inline void SampleSet2D::setTheta(int index, float angle)
{
    thetaAngles_[index] = angle;
    equalIntervalTheta_ = array_util::isEqualInterval(thetaAngles_);
}

inline void SampleSet2D::setPhi(int index, float angle)
{
    phiAngles_[index] = angle;
    equalIntervalPhi_ = array_util::isEqualInterval(phiAngles_);
}

inline Arrayf& SampleSet2D::getThetaArray() { return thetaAngles_; }
inline Arrayf& SampleSet2D::getPhiArray()   { return phiAngles_; }

inline const Arrayf& SampleSet2D::getThetaArray() const { return thetaAngles_; }
inline const Arrayf& SampleSet2D::getPhiArray()   const { return phiAngles_; }

inline int SampleSet2D::getNumTheta() const { return static_cast<int>(thetaAngles_.size()); }
inline int SampleSet2D::getNumPhi()   const { return static_cast<int>(phiAngles_.size()); }

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

inline int SampleSet2D::getNumWavelengths() const { return static_cast<int>(wavelengths_.size()); }

inline SourceType SampleSet2D::getSourceType() const { return sourceType_; }

inline void SampleSet2D::setSourceType(SourceType type) { sourceType_ = type; }

inline bool SampleSet2D::isIsotropic() const { return (phiAngles_.size() == 1); }

} // namespace lb

#endif // LIBBSDF_SAMPLE_SET_2D_H
