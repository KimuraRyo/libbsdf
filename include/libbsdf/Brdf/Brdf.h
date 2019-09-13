// =================================================================== //
// Copyright (C) 2014-2019 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#ifndef LIBBSDF_BRDF_H
#define LIBBSDF_BRDF_H

#include <libbsdf/Brdf/Sampler.h>
#include <libbsdf/Brdf/SampleSet.h>

#include <libbsdf/Common/Global.h>
#include <libbsdf/Common/Utility.h>
#include <libbsdf/Common/Log.h>

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

    /*! Gets the sample set including angles, wavelengths, and spectra. */
    SampleSet* getSampleSet();

    /*! Gets the sample set including angles, wavelengths, and spectra. */
    const SampleSet* getSampleSet() const;

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
     * \param verbose If this parameter is true, all messages of lb::Log::Level::WARN_MSG are output.
     */
    virtual bool validate(bool verbose = false) const = 0;

    /*!
     * Expands minimum angles to 0 and maximum angles to MAX_ANGLE,
     * and constructs the extrapolated sample set.
     */
    virtual bool expandAngles(bool angle0Expanded = true,
                              bool angle1Expanded = true,
                              bool angle2Expanded = true,
                              bool angle3Expanded = true) = 0;

    /*! Clamps all angles to minimum and maximum values of each coordinate system. */
    virtual void clampAngles() = 0;

    /*!
     * \brief Initializes all spectra of a BRDF using another BRDF.
     *
     * BRDFs must have the same wavelengths.
     */
    template <typename InterpolatorT>
    static bool initializeSpectra(const Brdf& baseBrdf, Brdf* brdf);

protected:
    /*! This attribute holds the sample set including angles, wavelengths, and spectra. */
    SampleSet* samples_;

    SourceType sourceType_; /*!< The data type of source. */

private:
    /*! Copy operator is disabled. */
    Brdf& operator=(const Brdf&);
};

inline       SampleSet* Brdf::getSampleSet()       { return samples_; }
inline const SampleSet* Brdf::getSampleSet() const { return samples_; }

inline SourceType Brdf::getSourceType() const { return sourceType_; }

inline void Brdf::setSourceType(SourceType type) { sourceType_ = type; }

template <typename InterpolatorT>
bool Brdf::initializeSpectra(const Brdf& baseBrdf, Brdf* brdf)
{
    lbTrace << "[Brdf::initializeSpectra]";

    const SampleSet* baseSs = baseBrdf.getSampleSet();
    SampleSet* ss = brdf->getSampleSet();

    if (!hasSameColor(*baseSs, *ss)) {
        lbError << "[Brdf::initializeSpectra] Color models or wavelengths do not match.";
        return false;
    }

    for (int i0 = 0; i0 < ss->getNumAngles0(); ++i0) {
    for (int i1 = 0; i1 < ss->getNumAngles1(); ++i1) {
    for (int i2 = 0; i2 < ss->getNumAngles2(); ++i2) {
    for (int i3 = 0; i3 < ss->getNumAngles3(); ++i3) {
        Vec3 inDir, outDir;
        brdf->getInOutDirection(i0, i1, i2, i3, &inDir, &outDir);

        Spectrum sp;
        Sampler::getSpectrum<InterpolatorT>(baseBrdf, inDir, outDir, &sp);

        ss->setSpectrum(i0, i1, i2, i3, sp);
    }}}}

    return true;
}

} // namespace lb

#endif // LIBBSDF_BRDF_H
