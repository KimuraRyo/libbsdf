// =================================================================== //
// Copyright (C) 2015-2019 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#ifndef LIBBSDF_RANDOM_SAMPLE_SET_H
#define LIBBSDF_RANDOM_SAMPLE_SET_H

#include <map>
#include <vector>

#include <libbsdf/Common/Utility.h>

namespace lb {

/*!
 * \class   RandomSampleSet
 * \brief   The RandomSampleSet class provides the BRDF using random sample points.
 *
 * The data structure consists of pairs of angles and a spectrum.
 * The coordinate system of angles is defined by \a CoordSysT.
 */
template <typename CoordSysT>
class RandomSampleSet
{
public:
    using AngleList = std::vector<float>;
    using SampleMap = std::map<AngleList, Spectrum, std::less<AngleList>,
                               Eigen::aligned_allocator<std::pair<const AngleList, Spectrum>>>;

    virtual ~RandomSampleSet() {}

    /*! Gets random sample points. */
    SampleMap& getSampleMap();

    /*! Gets random sample points. */
    const SampleMap& getSampleMap() const;

    /*! Finds the nearest sample and returns the spectrum. */
    const Spectrum& findSpectrumOfNearestSample(const AngleList& angles,
                                                bool             reciprocity = false) const;

    /*! Estimates the spectrum of a set of angles using a coordinate system (\a LocalCoordSysT). */
    template <typename LocalCoordSysT>
    const Spectrum& estimateSpectrum(const AngleList&   angles,
                                     float              weight0 = 1.0f,
                                     float              weight1 = 1.0f,
                                     float              weight2 = 1.0f,
                                     float              weight3 = 1.0f) const;

protected:
    SampleMap sampleMap_; /*!< Random sample points. */
};

template <typename CoordSysT>
inline typename RandomSampleSet<CoordSysT>::SampleMap& RandomSampleSet<CoordSysT>::getSampleMap()
{
    return sampleMap_;
}

template <typename CoordSysT>
inline const typename RandomSampleSet<CoordSysT>::SampleMap& RandomSampleSet<CoordSysT>::getSampleMap() const
{
    return sampleMap_;
}

template <typename CoordSysT>
const Spectrum& RandomSampleSet<CoordSysT>::findSpectrumOfNearestSample(const AngleList& angles,
                                                                        bool             reciprocity) const
{
    using std::acos;

    const Spectrum* sp = 0;

    Vec3 inDir, outDir;
    CoordSysT::toXyz(angles.at(0), angles.at(1), angles.at(2), angles.at(3),
                     &inDir, &outDir);

    Vec3::Scalar minAngleDiff = std::numeric_limits<Vec3::Scalar>::max();

    for (auto it = sampleMap_.begin(); it != sampleMap_.end(); ++it) {
        const AngleList& sampleAngles = it->first;
        Vec3 sampleInDir, sampleOutDir;
        CoordSysT::toXyz(sampleAngles.at(0), sampleAngles.at(1), sampleAngles.at(2), sampleAngles.at(3),
                         &sampleInDir, &sampleOutDir);

        Vec3::Scalar inAngle = acos(inDir.dot(sampleInDir));
        Vec3::Scalar outAngle = acos(outDir.dot(sampleOutDir));
        Vec3::Scalar angleDiff = inAngle + outAngle;

        if (reciprocity) {
            Vec3::Scalar inOutAngle = acos(inDir.dot(sampleOutDir));
            Vec3::Scalar outInAngle = acos(outDir.dot(sampleInDir));
            Vec3::Scalar angleUsingReciprocity = inOutAngle + outInAngle;

            if (angleUsingReciprocity < angleDiff) {
                angleDiff = angleUsingReciprocity;
            }
        }

        if (angleDiff == 0) {
            return it->second;
        }
        else if (angleDiff < minAngleDiff) {
            minAngleDiff = angleDiff;
            sp = &it->second;
        }
    }

    return *sp;
}

template <typename CoordSysT>
template <typename LocalCoordSysT>
const Spectrum& RandomSampleSet<CoordSysT>::estimateSpectrum(const AngleList&   angles,
                                                             float              weight0,
                                                             float              weight1,
                                                             float              weight2,
                                                             float              weight3) const
{
    using std::abs;

    float angle0, angle1, angle2, angle3;
    convertCoordinateSystem<CoordSysT, LocalCoordSysT>(
        angles.at(0), angles.at(1), angles.at(2), angles.at(3),
        &angle0, &angle1, &angle2, &angle3);

    const Spectrum* sp = 0;
    float minAngleDiff = std::numeric_limits<float>::max();

    for (auto it = sampleMap_.begin(); it != sampleMap_.end(); ++it) {
        const AngleList& sampleAngles = it->first;
        float sampleAngle0, sampleAngle1, sampleAngle2, sampleAngle3;
        convertCoordinateSystem<CoordSysT, LocalCoordSysT>(
            sampleAngles.at(0), sampleAngles.at(1), sampleAngles.at(2), sampleAngles.at(3),
            &sampleAngle0, &sampleAngle1, &sampleAngle2, &sampleAngle3);

        float angle0Diff = abs(angle0 - sampleAngle0);

        float angle1Diff = abs(angle1 - sampleAngle1);
        if (angle1Diff > PI_F) {
            angle1Diff = TAU_F - angle1Diff;
        }

        float angle2Diff = abs(angle2 - sampleAngle2);

        float angle3Diff = abs(angle3 - sampleAngle3);
        if (angle3Diff > PI_F) {
            angle3Diff = TAU_F - angle3Diff;
        }

        float angleDiff = weight0 * angle0Diff
                        + weight1 * angle1Diff
                        + weight2 * angle2Diff
                        + weight3 * angle3Diff;

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
