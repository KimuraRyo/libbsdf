// =================================================================== //
// Copyright (C) 2014-2015 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#ifndef LIBBSDF_HALF_DIFFERENCE_COORDINATE_SYSTEM_H
#define LIBBSDF_HALF_DIFFERENCE_COORDINATE_SYSTEM_H

#include <Eigen/Geometry>

#include <libbsdf/Common/SphericalCoordinateSystem.h>
#include <libbsdf/Common/Utility.h>

namespace lb {

/*!
 * \struct  HalfDifferenceCoordinateSystem
 * \brief   The HalfDifferenceCoordinateSystem struct provides the functions of a half difference coordinate system.
 *
 * See Rusinkiewicz, S. 1998. A New Change of Variables for Efficient BRDF Representation.
 *
 * The coordinate system has four angle parameters.
 *   - \a halfTheta: the polar angle of a halfway vector
 *   - \a halfPhi: the azimuthal angle of a halfway vector
 *   - \a diffTheta: the polar angle of a difference vector
 *   - \a diffPhi: the azimuthal angle of a difference vector
 *
 * \a diff is an abbreviation for difference. \a halfPhi is not used for isotropic BRDFs.
 */
struct HalfDifferenceCoordinateSystem
{
    /*!
     * Converts from four angles to incoming and outgoing directions and
     * assigns them to \a inDir and \a outDir.
     */
    static void toXyz(float halfTheta, float halfPhi,
                      float diffTheta, float diffPhi,
                      Vec3* inDir, Vec3* outDir);

    /*!
     * Converts from incoming and outgoing directions to four angles and
     * assigns them to \a halfTheta, \a halfPhi, \a diffTheta, and \a diffPhi.
     */
    static void fromXyz(const Vec3& inDir, const Vec3& outDir,
                        float* halfTheta, float* halfPhi,
                        float* diffTheta, float* diffPhi);

    /*!
     * Converts from incoming and outgoing directions to three angles for isotropic data and
     * assigns them to \a halfTheta, \a diffTheta, and \a diffPhi.
     */
    static void fromXyz(const Vec3& inDir, const Vec3& outDir,
                        float* halfTheta,
                        float* diffTheta, float* diffPhi);

    static const std::string ANGLE0_NAME; /*!< This attribute holds the name of halfTheta. */
    static const std::string ANGLE1_NAME; /*!< This attribute holds the name of halfPhi. */
    static const std::string ANGLE2_NAME; /*!< This attribute holds the name of diffTheta. */
    static const std::string ANGLE3_NAME; /*!< This attribute holds the name of diffPhi. */

    static const float MAX_ANGLE0; /*!< This attribute holds the maximum value of halfTheta. */
    static const float MAX_ANGLE1; /*!< This attribute holds the maximum value of halfPhi. */
    static const float MAX_ANGLE2; /*!< This attribute holds the maximum value of diffTheta. */
    static const float MAX_ANGLE3; /*!< This attribute holds the maximum value of diffPhi. */
};

inline void HalfDifferenceCoordinateSystem::toXyz(float halfTheta, float halfPhi,
                                                  float diffTheta, float diffPhi,
                                                  Vec3* inDir, Vec3* outDir)
{
    Vec3 halfDir = SphericalCoordinateSystem::toXyz(halfTheta, halfPhi);
    Vec3 diffDir = SphericalCoordinateSystem::toXyz(diffTheta, diffPhi);

    Vec2f rotThVec = Eigen::Rotation2D<Vec2f::Scalar>(-halfTheta) * Vec2f(diffDir[0], diffDir[2]);
    Vec2f rotPhVec = Eigen::Rotation2D<Vec2f::Scalar>(-halfPhi) * Vec2f(rotThVec[0], diffDir[1]);
    *inDir = Vec3(rotPhVec[0], rotPhVec[1], rotThVec[1]);

    *outDir = reflect(*inDir, halfDir);
}

inline void HalfDifferenceCoordinateSystem::fromXyz(const Vec3& inDir, const Vec3& outDir,
                                                    float* halfTheta, float* halfPhi,
                                                    float* diffTheta, float* diffPhi)
{
    Vec3 halfDir = (inDir + outDir).normalized();
    SphericalCoordinateSystem::fromXyz(halfDir, halfTheta, halfPhi);

    Vec2f rotPhVec = Eigen::Rotation2D<Vec2f::Scalar>(*halfPhi) * Vec2f(inDir[0], inDir[1]);
    Vec2f rotThVec = Eigen::Rotation2D<Vec2f::Scalar>(*halfTheta) * Vec2f(rotPhVec[0], inDir[2]);
    Vec3 diffDir(rotThVec[0], rotPhVec[1], rotThVec[1]);
    diffDir.normalize();
    SphericalCoordinateSystem::fromXyz(diffDir, diffTheta, diffPhi);
}

inline void HalfDifferenceCoordinateSystem::fromXyz(const Vec3& inDir, const Vec3& outDir,
                                                    float* halfTheta,
                                                    float* diffTheta, float* diffPhi)
{
    Vec3 halfDir = (inDir + outDir).normalized();
    float halfPhi; // halfPhi is 0 for isotorpic data.
    SphericalCoordinateSystem::fromXyz(halfDir, halfTheta, &halfPhi);

    Vec2f rotPhVec = Eigen::Rotation2D<Vec2f::Scalar>(-halfPhi) * Vec2f(inDir[0], inDir[1]);
    Vec2f rotThVec = Eigen::Rotation2D<Vec2f::Scalar>(*halfTheta) * Vec2f(rotPhVec[0], inDir[2]);
    Vec3 diffDir(rotThVec[0], rotPhVec[1], rotThVec[1]);
    diffDir.normalize();
    SphericalCoordinateSystem::fromXyz(diffDir, diffTheta, diffPhi);
}

} // namespace lb

#endif // LIBBSDF_HALF_DIFFERENCE_COORDINATE_SYSTEM_H
