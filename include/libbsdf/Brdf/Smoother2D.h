// =================================================================== //
// Copyright (C) 2017-2020 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#ifndef LIBBSDF_SMOOTHER2D_H
#define LIBBSDF_SMOOTHER2D_H

#include <set>

#include <libbsdf/Brdf/SampleSet2D.h>

namespace lb {

/*!
 * \class   Smoother2D
 * \brief   The Smoother2D class provides functions to smooth a 2D sample array by inserting angles.
 */
class Smoother2D
{
public:
    /*! Constructs the Smoother2D for a 2D sample array. */
    explicit Smoother2D(SampleSet2D* samples);

    /*! Smooths a 2D sample array by increasing sample points. */
    void smooth();

    /*!
     * Gets the threshold to insert sample points.
     * This is the difference between linear and smooth interpolation.
     */
    float getDiffThreshold() const;

    /*!
     * Sets the threshold to insert sample points.
     * This is the difference between linear and smooth interpolation.
     */
    void setDiffThreshold(float threshold);

    /*! Gets the maximum number of iterations of the division of angle0. */
    int getMaxIteration0() const;

    /*! Sets the maximum number of iterations of the division of angle0. */
    void setMaxIteration0(int maxIteration);

    /*! Gets the maximum number of iterations of the division of angle1. */
    int getMaxIteration1() const;

    /*! Sets the maximum number of iterations of the division of angle1. */
    void setMaxIteration1(int maxIteration);

private:
    void initializeAngles();

    bool insertAngle0();
    bool insertAngle1();

    bool insertAngle(std::set<Arrayf::Scalar>&  angleSet,
                     int                        angleSuffix,
                     const Vec2f&               angles,
                     const Vec2f&               nextAngles);

    void updateSamples();

    SampleSet2D* samples_;

    /*!
     * The difference between linear and smooth interpolation.
     * This attribute holds whether an angle is inserted.
     */
    float   diffThreshold_;

    int     maxIteration0_; /*!< The maximum number of iterations of the division of angle0. */
    int     maxIteration1_; /*!< The maximum number of iterations of the division of angle1. */

    float   minAngleInterval_; /*!< The minimum interval of divided angles. */

    std::set<Arrayf::Scalar> angles0_; /*!< The angle array to insert sample points. */
    std::set<Arrayf::Scalar> angles1_; /*!< The angle array to insert sample points. */
};

inline float Smoother2D::getDiffThreshold() const { return diffThreshold_; }
inline void Smoother2D::setDiffThreshold(float threshold) { diffThreshold_ = threshold; }

inline int Smoother2D::getMaxIteration0() const { return maxIteration0_; }
inline int Smoother2D::getMaxIteration1() const { return maxIteration1_; }

inline void Smoother2D::setMaxIteration0(int maxIteration) { maxIteration0_ = maxIteration; }
inline void Smoother2D::setMaxIteration1(int maxIteration) { maxIteration1_ = maxIteration; }

} // namespace lb

#endif // LIBBSDF_SMOOTHER2D_H
