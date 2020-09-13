// =================================================================== //
// Copyright (C) 2020 Kimura Ryo                                       //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#ifndef LIBBSDF_OPTIMIZER_H
#define LIBBSDF_OPTIMIZER_H

#include <set>

#include <libbsdf/Brdf/Brdf.h>

namespace lb {

/*!
 * \class   Optimizer
 * \brief   The Optimizer class provides functions to optimize a BRDF by removing angles.
 */
class Optimizer
{
public:
    /*! Constructs the optimizer for BRDF. */
    explicit Optimizer(Brdf* brdf,
                       float diffThreshold = 0.001f,
                       float ratioThreshold = 0.001f);

    /*! Optimizes a BRDF by removing sample points. */
    void optimize();

    /*!
     * Gets the difference threshold to judge sample points.
     * Metrics is the difference between a sampling point and linearly interpolated one.
     */
    float getDiffThreshold() const;

    /*!
     * Sets the difference threshold to judge sample points.
     * Metrics is the difference between a sampling point and linearly interpolated one.
     */
    void setDiffThreshold(float threshold);

    /*!
     * Gets the ratio threshold to judge sample points.
     * Metrics is the ratio between a sampling point and variation by linear interpolation.
     */
    float getRatioThreshold() const;

    /*!
     * Sets the ratio threshold to judge sample points.
     * Metrics is the ratio between a sampling point and variation by linear interpolation.
     */
    void setRatioThreshold(float threshold);

private:
    void setupAngles0();
    void setupAngles1();
    void setupAngles2();
    void setupAngles3();

    bool isExtraAngle(float         angle,
                      float         prevAngle,
                      float         nextAngle,
                      bool          prevDownward,
                      bool          nextDownward,
                      const Arrayd& sp,
                      const Arrayd& prevSp,
                      const Arrayd& nextSp);

    void updateBrdf();

    Brdf* brdf_;

    /*!
     * The difference between a sampling point and linearly interpolated one.
     * This attribute holds whether an angle is removed.
     */
    float diffThreshold_;

    /*!
     * The ratio between a sampling point and variation by linear interpolation.
     * This attribute holds whether an angle is removed.
     */
    float ratioThreshold_;

    std::set<Arrayf::Scalar> angles0_; /*!< The angle array to construct the new sample set. */
    std::set<Arrayf::Scalar> angles1_; /*!< The angle array to construct the new sample set. */
    std::set<Arrayf::Scalar> angles2_; /*!< The angle array to construct the new sample set. */
    std::set<Arrayf::Scalar> angles3_; /*!< The angle array to construct the new sample set. */
};

inline float Optimizer::getDiffThreshold() const { return diffThreshold_; }
inline void Optimizer::setDiffThreshold(float threshold) { diffThreshold_ = threshold; }

inline float Optimizer::getRatioThreshold() const { return ratioThreshold_; }
inline void Optimizer::setRatioThreshold(float threshold) { ratioThreshold_ = threshold; }

} // namespace lb

#endif // LIBBSDF_OPTIMIZER_H
