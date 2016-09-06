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
#include <libbsdf/Brdf/SpecularCoordinatesBrdf.h>
#include <libbsdf/Common/Global.h>

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
    typedef std::vector<float> AngleList;
    typedef std::map<AngleList, Spectrum, std::less<AngleList>,
                     Eigen::aligned_allocator<std::pair<AngleList, Spectrum> > > SampleMap;

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

    /*! Sets up a BRDF using random sample points. */
    void setupBrdf(SphericalCoordinatesBrdf* brdf);

    /*! Sets up a BRDF using random sample points. */
    void setupBrdf(SpecularCoordinatesBrdf* brdf,
                   float                    weight0 = 1.0f,
                   float                    weight1 = 1.0f,
                   float                    weight2 = 1.0f,
                   float                    weight3 = 1.0f);

private:
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
                                                                        bool              reciprocity) const
{
    using std::acos;

    const Spectrum* sp;

    Vec3 inDir, outDir;
    CoordSysT::toXyz(angles.at(0), angles.at(1), angles.at(2), angles.at(3),
                     &inDir, &outDir);

    float minAngleDiff = std::numeric_limits<float>::max();

    for (auto it = sampleMap_.begin(); it != sampleMap_.end(); ++it) {
        const AngleList& sampleAngles = it->first;
        Vec3 sampleInDir, sampleOutDir;
        CoordSysT::toXyz(sampleAngles.at(0), sampleAngles.at(1), sampleAngles.at(2), sampleAngles.at(3),
                         &sampleInDir, &sampleOutDir);

        float inAngle = acos(inDir.dot(sampleInDir));
        float outAngle = acos(outDir.dot(sampleOutDir));
        float angleDiff = inAngle + outAngle;

        if (reciprocity) {
            float inOutAngle = acos(inDir.dot(sampleOutDir));
            float outInAngle = acos(outDir.dot(sampleInDir));
            float angleUsingReciprocity = inOutAngle + outInAngle;

            if (angleUsingReciprocity < angleDiff) {
                angleDiff = angleUsingReciprocity;
            }
        }

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

    const Spectrum* sp;
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
            angle1Diff = 2.0f * PI_F - angle1Diff;
        }

        float angle2Diff = abs(angle2 - sampleAngle2);

        float angle3Diff = abs(angle3 - sampleAngle3);
        if (angle3Diff > PI_F) {
            angle3Diff = 2.0f * PI_F - angle3Diff;
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

template <typename CoordSysT>
void RandomSampleSet<CoordSysT>::setupBrdf(SphericalCoordinatesBrdf* brdf)
{
    for (int inThIndex = 0;  inThIndex < brdf->getNumInTheta();   ++inThIndex)  {
    for (int inPhIndex = 0;  inPhIndex < brdf->getNumInPhi();     ++inPhIndex)  {
    for (int outThIndex = 0; outThIndex < brdf->getNumOutTheta(); ++outThIndex) {
        RandomSampleSet::AngleList angles;
        SampleMap::iterator it;
        #pragma omp parallel for private(angles, it)
        for (int outPhIndex = 0; outPhIndex < brdf->getNumOutPhi(); ++outPhIndex) {
            angles.resize(4);
            angles.at(0) = brdf->getInTheta(inThIndex);
            angles.at(1) = brdf->getInPhi(inPhIndex);
            angles.at(2) = brdf->getOutTheta(outThIndex);
            angles.at(3) = brdf->getOutPhi(outPhIndex);

            convertCoordinateSystem<SphericalCoordinateSystem, CoordSysT>(
                angles.at(0), angles.at(1), angles.at(2), angles.at(3),
                &angles[0], &angles[1], &angles[2], &angles[3]);

            it = sampleMap_.find(angles);
            if (it != sampleMap_.end()) {
                brdf->setSpectrum(inThIndex, inPhIndex, outThIndex, outPhIndex,
                                  it->second);
            }
            else {
                brdf->setSpectrum(inThIndex, inPhIndex, outThIndex, outPhIndex,
                                  findSpectrumOfNearestSample(angles, false));
                //brdf->setSpectrum(inThIndex, inPhIndex, outThIndex, outPhIndex,
                //                  estimateSpectrum<CoordSysT>(angles));
            }
        }
    }}}

    brdf->getSampleSet()->updateAngleAttributes();
}

template <typename CoordSysT>
void RandomSampleSet<CoordSysT>::setupBrdf(SpecularCoordinatesBrdf* brdf,
                                           float                    weight0,
                                           float                    weight1,
                                           float                    weight2,
                                           float                    weight3)
{
    for (int inThIndex = 0; inThIndex < brdf->getNumInTheta();   ++inThIndex) {
    for (int inPhIndex = 0; inPhIndex < brdf->getNumInPhi();     ++inPhIndex) {
    for (int spThIndex = 0; spThIndex < brdf->getNumSpecTheta(); ++spThIndex) {
        RandomSampleSet::AngleList angles;
        SampleMap::iterator it;
        float w3;
        #pragma omp parallel for private(angles, it, w3)
        for (int spPhIndex = 0; spPhIndex < brdf->getNumSpecPhi(); ++spPhIndex) {
            angles.resize(4);
            angles.at(0) = brdf->getInTheta(inThIndex);
            angles.at(1) = brdf->getInPhi(inPhIndex);
            angles.at(2) = brdf->getSpecTheta(spThIndex);
            angles.at(3) = brdf->getSpecPhi(spPhIndex);

            convertCoordinateSystem<SpecularCoordinateSystem, CoordSysT>(
                angles.at(0), angles.at(1), angles.at(2), angles.at(3),
                &angles[0], &angles[1], &angles[2], &angles[3]);

            it = sampleMap_.find(angles);
            if (it != sampleMap_.end()) {
                brdf->setSpectrum(inThIndex, inPhIndex, spThIndex, spPhIndex,
                                  it->second);
            }
            else {
                // Modify a weight coefficient.
                w3 = weight3 * hermiteInterpolation3(1.0f, 1.0f / weight3, brdf->getSpecTheta(spThIndex) / PI_2_F);
                Spectrum sp = estimateSpectrum<SpecularCoordinateSystem>(angles, weight0, weight1, weight2, w3);
                brdf->setSpectrum(inThIndex, inPhIndex, spThIndex, spPhIndex, sp);
            }
        }
    }}}

    brdf->getSampleSet()->updateAngleAttributes();
}

} // namespace lb

#endif // LIBBSDF_RANDOM_SAMPLE_SET_H
