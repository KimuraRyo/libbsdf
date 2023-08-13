// =================================================================== //
// Copyright (C) 2017-2023 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#ifndef LIBBSDF_SMOOTHER_H
#define LIBBSDF_SMOOTHER_H

#include <set>

#include <libbsdf/Brdf/Brdf.h>

namespace lb {

/*!
 * \class   Smoother
 * \brief   The Smoother class provides functions to smooth a BRDF by inserting angles.
 */
class Smoother
{
public:
    /*! Constructs the smoother for BRDF. */
    explicit Smoother(Brdf* brdf);

    /*! Smooths a BRDF by increasing sample points. */
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

    /*! Gets the maximum number of iterations of the division of angle2. */
    int getMaxIteration2() const;

    /*! Sets the maximum number of iterations of the division of angle2. */
    void setMaxIteration2(int maxIteration);

    /*! Gets the maximum number of iterations of the division of angle3. */
    int getMaxIteration3() const;

    /*! Sets the maximum number of iterations of the division of angle3. */
    void setMaxIteration3(int maxIteration);

    /*! Gets the maximum specular polar angle to avoid smoothing. */
    double getSpecularPolarRegion() const;

    /*! Sets the maximum specular polar angle to avoid smoothing. */
    void setSpecularPolarRegion(double angle);

private:
    void initializeAngles();

    bool insertAngle0();
    bool insertAngle1();
    bool insertAngle2();
    bool insertAngle3();

    bool insertAngle(std::set<double>& angleSet,
                     int               angleSuffix,
                     const Vec4&       angles,
                     const Vec4&       nextAngles);

    void updateBrdf();

    Brdf* brdf_;

    /*!
     * The difference between linear and smooth interpolation.
     * This attribute holds whether an angle is inserted.
     */
    float diffThreshold_;

    int maxIteration0_; /*!< The maximum number of iterations of the division of angle0. */
    int maxIteration1_; /*!< The maximum number of iterations of the division of angle1. */
    int maxIteration2_; /*!< The maximum number of iterations of the division of angle2. */
    int maxIteration3_; /*!< The maximum number of iterations of the division of angle3. */

    double minAngleInterval_; /*!< The minimum interval of divided angles. */

    /*!
     * The maximum specular polar angle to avoid smoothing.
     * If this attribute is 0, all regions are smoothed.
     * \a brdf_ must be \a SpecularCoordinatesBrdf.
     */
    double specularPolarRegion_;

    std::set<double> angles0_; /*!< The angle array to insert sample points. */
    std::set<double> angles1_; /*!< The angle array to insert sample points. */
    std::set<double> angles2_; /*!< The angle array to insert sample points. */
    std::set<double> angles3_; /*!< The angle array to insert sample points. */
};

inline float Smoother::getDiffThreshold() const { return diffThreshold_; }
inline void Smoother::setDiffThreshold(float threshold) { diffThreshold_ = threshold; }

inline int Smoother::getMaxIteration0() const { return maxIteration0_; }
inline int Smoother::getMaxIteration1() const { return maxIteration1_; }
inline int Smoother::getMaxIteration2() const { return maxIteration2_; }
inline int Smoother::getMaxIteration3() const { return maxIteration3_; }

inline void Smoother::setMaxIteration0(int maxIteration) { maxIteration0_ = maxIteration; }
inline void Smoother::setMaxIteration1(int maxIteration) { maxIteration1_ = maxIteration; }
inline void Smoother::setMaxIteration2(int maxIteration) { maxIteration2_ = maxIteration; }
inline void Smoother::setMaxIteration3(int maxIteration) { maxIteration3_ = maxIteration; }

inline double Smoother::getSpecularPolarRegion() const { return specularPolarRegion_; }
inline void   Smoother::setSpecularPolarRegion(double angle) { specularPolarRegion_ = angle; }

} // namespace lb

#endif // LIBBSDF_SMOOTHER_H
