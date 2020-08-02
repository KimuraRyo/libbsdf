// =================================================================== //
// Copyright (C) 2014-2020 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#ifndef LIBBSDF_BRDF_H
#define LIBBSDF_BRDF_H

#include <string>

#include <libbsdf/Brdf/SampleSet.h>

namespace lb {

/*!
 * \class   Brdf
 * \brief   The Brdf class provides the BRDF data and sampling functions.
 *
 * The BRDF data consists of angles, wavelengths, spectra, and coordinate system.
 * The data structure is defined in lb::SampleSet. The functions depending on coordinate
 * system are implemented in lb::Brdf and derived classes.
 */
class Brdf
{
public:
    /*! Constructs a BRDF. */
    Brdf(int        numAngles0,
         int        numAngles1,
         int        numAngles2,
         int        numAngles3,
         ColorModel colorModel = RGB_MODEL,
         int        numWavelengths = 3);

    /*! Constructs an empty BRDF. Brdf::samples_ must be initialized in a derived class. */
    Brdf();

    /*! Copies and constructs a BRDF. */
    Brdf(const Brdf& brdf);

    virtual ~Brdf();

    /*! Virtual copy constructor. */
    virtual Brdf* clone() const = 0;

    /*! Initializes spectra using lb::Brdf. */
    void initializeSpectra(const Brdf& brdf);

    /*! Gets the sample set including angles, wavelengths, and spectra. */
    SampleSet* getSampleSet();

    /*! Gets the sample set including angles, wavelengths, and spectra. */
    const SampleSet* getSampleSet() const;

    /*! Gets the source type. */
    ReductionType getReductionType() const;

    /*! Sets the source type. */
    void setReductionType(ReductionType type);

    /*! Gets the source type. */
    SourceType getSourceType() const;

    /*! Sets the source type. */
    void setSourceType(SourceType type);

    /*! Gets the spectrum of the BRDF at incoming and outgoing directions. */
    virtual Spectrum getSpectrum(const Vec3& inDir, const Vec3& outDir) const = 0;

    /*! Gets the value of the BRDF at incoming and outgoing directions and the index of wavelength. */
    virtual float getValue(const Vec3& inDir, const Vec3& outDir, int wavelengthIndex) const = 0;

    /*!
     * Computes incoming and outgoing directions of a Cartesian coordinate system
     * using a set of angle indices.
     */
    virtual void getInOutDirection(int      index0,
                                   int      index1,
                                   int      index2,
                                   int      index3,
                                   Vec3*    inDir,
                                   Vec3*    outDir) const = 0;

    /*!
     * Converts from four angles to incoming and outgoing directions and
     * assigns them to \a inDir and \a outDir.
     */
    virtual void toXyz(float angle0,
                       float angle1,
                       float angle2,
                       float angle3,
                       Vec3* inDir,
                       Vec3* outDir) const = 0;

    /*!
     * Converts from incoming and outgoing directions to four angles and
     * assigns them to \a angle0, \a angle1, \a angle2, and \a angle3.
     */
    virtual void fromXyz(const Vec3&    inDir,
                         const Vec3&    outDir,
                         float*         angle0,
                         float*         angle1,
                         float*         angle2,
                         float*         angle3) const = 0;

    /*!
     * Converts from incoming and outgoing directions to three angles for an isotropic BRDF and
     * assigns them to \a angle0, \a angle2, and \a angle3.
     */
    virtual void fromXyz(const Vec3&    inDir,
                         const Vec3&    outDir,
                         float*         angle0,
                         float*         angle2,
                         float*         angle3) const = 0;

    virtual std::string getAngle0Name() const = 0; /*!< Gets the name of angle0. */
    virtual std::string getAngle1Name() const = 0; /*!< Gets the name of angle1. */
    virtual std::string getAngle2Name() const = 0; /*!< Gets the name of angle2. */
    virtual std::string getAngle3Name() const = 0; /*!< Gets the name of angle3. */

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
    virtual bool validate(bool verbose = false) const = 0;

    /*!
     * Expands minimum angles to MIN_ANGLE and maximum angles to MAX_ANGLE,
     * and constructs the extrapolated sample set.
     */
    virtual bool expandAngles(bool angle0Expanded = true,
                              bool angle1Expanded = true,
                              bool angle2Expanded = true,
                              bool angle3Expanded = true) = 0;

    /*! Clamps all angles to minimum and maximum values of each coordinate system. */
    virtual void clampAngles() = 0;

    /*! Gets the name of BRDF. */
    std::string& getName();

    /*! Gets the name of BRDF. */
    const std::string& getName() const;

    /*! Sets the name of BRDF. */
    void setName(const std::string& name);

protected:
    /*! This attribute holds the sample set including angles, wavelengths, and spectra. */
    SampleSet* samples_;

    ReductionType   reductionType_; /*!< The reduction type of data. */
    SourceType      sourceType_;    /*!< The data type of source. */

private:
    /*! Copy operator is disabled. */
    Brdf& operator=(const Brdf&);

    std::string name_; /*! The name of BRDF. */
};

inline       SampleSet* Brdf::getSampleSet()       { return samples_; }
inline const SampleSet* Brdf::getSampleSet() const { return samples_; }

inline ReductionType Brdf::getReductionType() const { return reductionType_; }

inline void Brdf::setReductionType(ReductionType type) { reductionType_ = type; }

inline SourceType Brdf::getSourceType() const { return sourceType_; }

inline void Brdf::setSourceType(SourceType type) { sourceType_ = type; }

inline       std::string& Brdf::getName()       { return name_; }
inline const std::string& Brdf::getName() const { return name_; }

} // namespace lb

#endif // LIBBSDF_BRDF_H
