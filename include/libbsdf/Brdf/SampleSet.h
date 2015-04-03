// =================================================================== //
// Copyright (C) 2014-2015 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#ifndef LIBBSDF_SAMPLE_SET_H
#define LIBBSDF_SAMPLE_SET_H

#include <cassert>

#include <libbsdf/Common/Spectrum.h>
#include <libbsdf/Common/Vector.h>

namespace lb {

/*!
 * \class   SampleSet
 * \brief   The SampleSet class provides the BRDF data structure.
 *
 * The data structure consists of arrays of angles, wavelengths, and spectra.
 * The arrays of angles and wavelengths must be sorted in ascending order.
 *
 * A sample point is defined with four angles.
 *   - angle0 (e.g. incoming polar angle of a spherical coordinate system)
 *   - angle1 (e.g. incoming azimuthal angle of a spherical coordinate system)
 *   - angle2 (e.g. outgoing polar angle of a spherical coordinate system)
 *   - angle3 (e.g. outgoing azimuthal angle of a spherical coordinate system)
 */
class SampleSet
{
public:
    /*! Constructs the sample points of a BRDF. */
    SampleSet(int               numAngles0,
              int               numAngles1,
              int               numAngles2,
              int               numAngles3,
              ColorModel::Type  colorModel = ColorModel::RGB,
              int               numWavelengths = 3);

    /*! Sets spectra and angles to zero. */
    void zeroOutSamples();

    /*! Gets the spectrum at a set of angle indices. */
    Spectrum& getSpectrum(int index0, int index1, int index2, int index3);

    /*! Gets the spectrum at a set of angle indices. */
    const Spectrum& getSpectrum(int index0, int index1, int index2, int index3) const;

    /*! Gets the spectrum at an index. */
    Spectrum& getSpectrum(int index);

    /*! Gets the spectrum at an index. */
    const Spectrum& getSpectrum(int index) const;

    /*! Sets the spectrum at a set of angle indices. */
    void setSpectrum(int index0, int index1, int index2, int index3,
                     const Spectrum& spectrum);

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

    float getAngle0(int index) const; /*!< Gets the angle0 at an index. */
    float getAngle1(int index) const; /*!< Gets the angle1 at an index. */
    float getAngle2(int index) const; /*!< Gets the angle2 at an index. */
    float getAngle3(int index) const; /*!< Gets the angle3 at an index. */

    void setAngle0(int index, float angle); /*!< Sets the angle0 at an index. */
    void setAngle1(int index, float angle); /*!< Sets the angle1 at an index. */
    void setAngle2(int index, float angle); /*!< Sets the angle2 at an index. */
    void setAngle3(int index, float angle); /*!< Sets the angle3 at an index. */

    Arrayf& getAngles0(); /*!< Gets the array of angles0. */
    Arrayf& getAngles1(); /*!< Gets the array of angles1. */
    Arrayf& getAngles2(); /*!< Gets the array of angles2. */
    Arrayf& getAngles3(); /*!< Gets the array of angles3. */

    const Arrayf& getAngles0() const; /*!< Gets the array of angles0. */
    const Arrayf& getAngles1() const; /*!< Gets the array of angles1. */
    const Arrayf& getAngles2() const; /*!< Gets the array of angles2. */
    const Arrayf& getAngles3() const; /*!< Gets the array of angles3. */

    int getNumAngles0() const; /*!< Gets the number of angles0. */
    int getNumAngles1() const; /*!< Gets the number of angles1. */
    int getNumAngles2() const; /*!< Gets the number of angles2. */
    int getNumAngles3() const; /*!< Gets the number of angles3. */

    bool isEqualIntervalAngles0() const; /*!< Returns true if angles0 are set at equal intervals. */
    bool isEqualIntervalAngles1() const; /*!< Returns true if angles1 are set at equal intervals. */
    bool isEqualIntervalAngles2() const; /*!< Returns true if angles2 are set at equal intervals. */
    bool isEqualIntervalAngles3() const; /*!< Returns true if angles3 are set at equal intervals. */

    /*! Gets the color model. */
    ColorModel::Type getColorModel() const;

    /*! Checks the attributes whether angles are set at equal intervals */
    void checkEqualIntervalAngles();

    /*! Resizes spectra and angles. */
    void resizeSamples(int numAngles0, int numAngles1, int numAngles2, int numAngles3);

private:
    /*! Resizes each spectrum in spectra_. */
    void initializeSpectralSamples(int numWavelengths);

    /*! Gets the index of spectra from a set of angle indices. */
    int getIndex(int index0, int index1, int index2, int index3) const;

    SpectrumList spectra_; /*!< The list of spectrum for each pair of incoming and outgoing directions. */

    Arrayf wavelengths_; /*!< The array of wavelengths. */

    Arrayf angles0_; /*!< The array of angles0. */
    Arrayf angles1_; /*!< The array of angles1. */
    Arrayf angles2_; /*!< The array of angles2. */
    Arrayf angles3_; /*!< The array of angles3. */

    int numAngles0_; /*!< The number of angles0. */
    int numAngles1_; /*!< The number of angles1. */
    int numAngles2_; /*!< The number of angles2. */
    int numAngles3_; /*!< The number of angles3. */

    bool equalIntervalAngles0_; /*!< This attribute holds whether angles0 are set at equal intervals. */
    bool equalIntervalAngles1_; /*!< This attribute holds whether angles1 are set at equal intervals. */
    bool equalIntervalAngles2_; /*!< This attribute holds whether angles2 are set at equal intervals. */
    bool equalIntervalAngles3_; /*!< This attribute holds whether angles3 are set at equal intervals. */

    ColorModel::Type colorModel_; /*!< The color model of spectra. */
};

inline Spectrum& SampleSet::getSpectrum(int index0, int index1, int index2, int index3)
{
    return spectra_.at(getIndex(index0, index1, index2, index3));
}

inline const Spectrum& SampleSet::getSpectrum(int index0, int index1, int index2, int index3) const
{
    return spectra_.at(getIndex(index0, index1, index2, index3));
}

inline       Spectrum& SampleSet::getSpectrum(int index)       { return spectra_.at(index); }
inline const Spectrum& SampleSet::getSpectrum(int index) const { return spectra_.at(index); }

inline void SampleSet::setSpectrum(int index0, int index1, int index2, int index3,
                                   const Spectrum& spectrum)
{
    spectra_.at(getIndex(index0, index1, index2, index3)) = spectrum;
}

inline float SampleSet::getWavelength(int index) const
{
    return wavelengths_[index];
}

inline void SampleSet::setWavelength(int index, float wavelength)
{
    wavelengths_[index] = wavelength;
}

inline       Arrayf& SampleSet::getWavelengths()       { return wavelengths_; }
inline const Arrayf& SampleSet::getWavelengths() const { return wavelengths_; }

inline int SampleSet::getNumWavelengths() const
{
    return wavelengths_.size();
}

inline float SampleSet::getAngle0(int index) const { return angles0_[index]; }
inline float SampleSet::getAngle1(int index) const { return angles1_[index]; }
inline float SampleSet::getAngle2(int index) const { return angles2_[index]; }
inline float SampleSet::getAngle3(int index) const { return angles3_[index]; }

inline void SampleSet::setAngle0(int index, float angle) { angles0_[index] = angle; }
inline void SampleSet::setAngle1(int index, float angle) { angles1_[index] = angle; }
inline void SampleSet::setAngle2(int index, float angle) { angles2_[index] = angle; }
inline void SampleSet::setAngle3(int index, float angle) { angles3_[index] = angle; }

inline Arrayf& SampleSet::getAngles0() { return angles0_; }
inline Arrayf& SampleSet::getAngles1() { return angles1_; }
inline Arrayf& SampleSet::getAngles2() { return angles2_; }
inline Arrayf& SampleSet::getAngles3() { return angles3_; }

inline const Arrayf& SampleSet::getAngles0() const { return angles0_; }
inline const Arrayf& SampleSet::getAngles1() const { return angles1_; }
inline const Arrayf& SampleSet::getAngles2() const { return angles2_; }
inline const Arrayf& SampleSet::getAngles3() const { return angles3_; }

inline int SampleSet::getNumAngles0() const { return numAngles0_; }
inline int SampleSet::getNumAngles1() const { return numAngles1_; }
inline int SampleSet::getNumAngles2() const { return numAngles2_; }
inline int SampleSet::getNumAngles3() const { return numAngles3_; }

inline bool SampleSet::isEqualIntervalAngles0() const { return equalIntervalAngles0_; }
inline bool SampleSet::isEqualIntervalAngles1() const { return equalIntervalAngles1_; }
inline bool SampleSet::isEqualIntervalAngles2() const { return equalIntervalAngles2_; }
inline bool SampleSet::isEqualIntervalAngles3() const { return equalIntervalAngles3_; }

inline ColorModel::Type SampleSet::getColorModel() const
{
    return colorModel_;
}

inline int SampleSet::getIndex(int index0, int index1, int index2, int index3) const
{
    assert(index0 >= 0 && index1 >= 0 && index2 >= 0 && index3 >= 0);
    assert(index0 < numAngles0_ && index1 < numAngles1_ && index2 < numAngles2_ && index3 < numAngles3_);

    int index = index0
              + numAngles0_ * index1
              + numAngles0_ * numAngles1_ * index2
              + numAngles0_ * numAngles1_ * numAngles2_ * index3;
    return index;
}

} // namespace lb

#endif // LIBBSDF_SAMPLE_SET_H
