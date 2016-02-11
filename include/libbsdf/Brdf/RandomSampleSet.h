// =================================================================== //
// Copyright (C) 2015-2016 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#ifndef LIBBSDF_RANDOM_SAMPLE_SET_H
#define LIBBSDF_RANDOM_SAMPLE_SET_H

#include <map>
#include <vector>

#include <libbsdf/Brdf/SphericalCoordinatesBrdf.h>
#include <libbsdf/Common/SpecularCoordinateSystem.h>
#include <libbsdf/Common/Global.h>

namespace lb {

/*!
 * \class   RandomSampleSet
 * \brief   The RandomSampleSet class provides the BRDF using random sample points.
 *
 * The data structure consists of pairs of angles and a spectrum.
 * Angles are a spherical coordinate system.
 */
class RandomSampleSet
{
public:
    typedef std::vector<float> AngleList;
    typedef std::map<AngleList, Spectrum, std::less<AngleList>,
                     Eigen::aligned_allocator<std::pair<AngleList, Spectrum> > > SampleMap;

    /*! Gets random sample points. */
    SampleMap& getSampleMap();

    /*! Gets random sample points. */
    const SampleMap& getSampleMap() const;

    /*! Finds the nearest sample and returns the spectrum. */
    const Spectrum& getSpectrumOfNearestSample(const AngleList& angles,
                                               bool             useReciprocity = false) const;

    /*! Estimates the spectrum of a set of angles using a coordinate system. */
    template <typename CoordSysT>
    const Spectrum& estimateSpectrum(const AngleList& angles) const;

    /*! Sets up a BRDF using random sample points. */
    void setupBrdf(SphericalCoordinatesBrdf* brdf);

private:
    SampleMap sampleMap_; /*!< Random sample points. */
};

inline       RandomSampleSet::SampleMap& RandomSampleSet::getSampleMap()       { return sampleMap_; }
inline const RandomSampleSet::SampleMap& RandomSampleSet::getSampleMap() const { return sampleMap_; }

template <typename CoordSysT>
const Spectrum& RandomSampleSet::estimateSpectrum(const AngleList& angles) const
{
    using std::abs;

    float inTheta, inPhi, specTheta, specPhi;
    convertCoordinateSystem<SphericalCoordinateSystem, SpecularCoordinateSystem> (
        angles.at(0), angles.at(1), angles.at(2), angles.at(3),
        &inTheta, &inPhi, &specTheta, &specPhi);

    const Spectrum* sp;
    float minAngleDiff = std::numeric_limits<float>::max();

    for (auto it = sampleMap_.begin(); it != sampleMap_.end(); ++it) {
        const AngleList& sampleAngles = it->first;
        float sampleInTheta, sampleInPhi, sampleSpecTheta, sampleSpecPhi;
        convertCoordinateSystem<SphericalCoordinateSystem, SpecularCoordinateSystem>(
            sampleAngles.at(0), sampleAngles.at(1), sampleAngles.at(2), sampleAngles.at(3),
            &sampleInTheta, &sampleInPhi, &sampleSpecTheta, &sampleSpecPhi);

        float inThetaDiff = abs(inTheta - sampleInTheta);

        float inPhiDiff = abs(inPhi - sampleInPhi);
        if (inPhiDiff > PI_F) {
            inPhiDiff = 2.0f * PI_F - inPhiDiff;
        }

        float specThetaDiff = abs(specTheta - sampleSpecTheta);

        float specPhiDiff = abs(specPhi - sampleSpecPhi);
        if (specPhiDiff > PI_F) {
            specPhiDiff = 2.0f * PI_F - specPhiDiff;
        }

        float angleDiff = inThetaDiff + inPhiDiff + specThetaDiff + specPhiDiff;

        if (angleDiff == 0.0f) {
            return it->second;
        }
        else if (angleDiff < minAngleDiff) {
            minAngleDiff = angleDiff;
            sp = &it->second;
        }
    }

    return *sp;
}

} // namespace lb

#endif // LIBBSDF_RANDOM_SAMPLE_SET_H
