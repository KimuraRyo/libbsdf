// =================================================================== //
// Copyright (C) 2014-2015 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#ifndef LIBBSDF_BRDF_H
#define LIBBSDF_BRDF_H

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
    virtual ~Brdf();

    /*! Gets the spectrum of the BRDF at incoming and outgoing directions. */
    virtual Spectrum getSpectrum(const Vec3& inDir, const Vec3& outDir) const = 0;

    /*! Gets the value of the BRDF at incoming and outgoing directions and the index of wavelength. */
    virtual float getValue(const Vec3& inDir, const Vec3& outDir, int wavelengthIndex) const = 0;

    /*!
     * Computes incoming and outgoing directions of a Cartesian coordinate system
     * using a set of angle indices.
     */
    virtual void getInOutDirection(int index0, int index1, int index2, int index3,
                                   Vec3* inDir, Vec3* outDir) const = 0;

    /*! Gets the sample set including angles, wavelengths, and spectra. */
    SampleSet* getSampleSet();

    /*! Gets the sample set including angles, wavelengths, and spectra. */
    const SampleSet* getSampleSet() const;

    virtual std::string getAngle0Name() const = 0; /*!< Gets the name of angle0. */
    virtual std::string getAngle1Name() const = 0; /*!< Gets the name of angle1. */
    virtual std::string getAngle2Name() const = 0; /*!< Gets the name of angle2. */
    virtual std::string getAngle3Name() const = 0; /*!< Gets the name of angle3. */

    /*!
     * Expands minimum angles to 0 and maximum angles to MAX_ANGLE,
     * and constructs the extrapolated sample set.
     */
    virtual bool expandAngles() = 0;

    /*! Clamps all angles to minimum and maximum values of each coordinate system. */
    virtual void clampAngles() = 0;

protected:
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

    /*! This attribute holds the sample set including angles, wavelengths, and spectra. */
    SampleSet* samples_;

private:
    /*! Copy operator is disabled. */
    Brdf& operator=(const Brdf&);
};

inline       SampleSet* Brdf::getSampleSet()       { return samples_; }
inline const SampleSet* Brdf::getSampleSet() const { return samples_; }

} // namespace lb

#endif // LIBBSDF_BRDF_H
