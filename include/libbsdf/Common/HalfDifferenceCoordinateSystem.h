// =================================================================== //
// Copyright (C) 2014-2019 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#ifndef LIBBSDF_HALF_DIFFERENCE_COORDINATE_SYSTEM_H
#define LIBBSDF_HALF_DIFFERENCE_COORDINATE_SYSTEM_H

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
    template <typename ScalarT>
    static void toXyz(ScalarT   halfTheta,
                      ScalarT   halfPhi,
                      ScalarT   diffTheta,
                      ScalarT   diffPhi,
                      Vec3*     inDir,
                      Vec3*     outDir);

    /*!
     * Converts from incoming and outgoing directions to four angles and
     * assigns them to \a halfTheta, \a halfPhi, \a diffTheta, and \a diffPhi.
     */
    template <typename ScalarT>
    static void fromXyz(const Vec3& inDir,
                        const Vec3& outDir,
                        ScalarT*    halfTheta,
                        ScalarT*    halfPhi,
                        ScalarT*    diffTheta,
                        ScalarT*    diffPhi);

    /*!
     * Converts from incoming and outgoing directions to three angles for isotropic data and
     * assigns them to \a halfTheta, \a diffTheta, and \a diffPhi.
     */
    template <typename ScalarT>
    static void fromXyz(const Vec3& inDir,
                        const Vec3& outDir,
                        ScalarT*    halfTheta,
                        ScalarT*    diffTheta,
                        ScalarT*    diffPhi);

    static const std::string ANGLE0_NAME; /*!< This attribute holds the name of halfTheta. */
    static const std::string ANGLE1_NAME; /*!< This attribute holds the name of halfPhi. */
    static const std::string ANGLE2_NAME; /*!< This attribute holds the name of diffTheta. */
    static const std::string ANGLE3_NAME; /*!< This attribute holds the name of diffPhi. */

    static const float MIN_ANGLE0; /*!< This attribute holds the minimum value of halfTheta. */
    static const float MIN_ANGLE1; /*!< This attribute holds the minimum value of halfPhi. */
    static const float MIN_ANGLE2; /*!< This attribute holds the minimum value of diffTheta. */
    static const float MIN_ANGLE3; /*!< This attribute holds the minimum value of diffPhi. */

    static const float MAX_ANGLE0; /*!< This attribute holds the maximum value of halfTheta. */
    static const float MAX_ANGLE1; /*!< This attribute holds the maximum value of halfPhi. */
    static const float MAX_ANGLE2; /*!< This attribute holds the maximum value of diffTheta. */
    static const float MAX_ANGLE3; /*!< This attribute holds the maximum value of diffPhi. */
};

template <typename ScalarT>
void HalfDifferenceCoordinateSystem::toXyz(ScalarT  halfTheta,
                                           ScalarT  halfPhi,
                                           ScalarT  diffTheta,
                                           ScalarT  diffPhi,
                                           Vec3*    inDir,
                                           Vec3*    outDir)
{
    Vec3 halfDir = SphericalCoordinateSystem::toXyz(halfTheta, halfPhi);
    Vec3 diffDir = SphericalCoordinateSystem::toXyz(diffTheta, diffPhi);

    Vec2 rotThVec = Eigen::Rotation2D<Vec2::Scalar>(-halfTheta) * Vec2(diffDir[0], diffDir[2]);
    Vec2 rotPhVec = Eigen::Rotation2D<Vec2::Scalar>(halfPhi) * Vec2(rotThVec[0], diffDir[1]);
    *inDir = Vec3(static_cast<Vec3::Scalar>(rotPhVec[0]),
                  static_cast<Vec3::Scalar>(rotPhVec[1]),
                  static_cast<Vec3::Scalar>(rotThVec[1]));

    (*inDir)[2] = std::max((*inDir)[2], Vec3::Scalar(0));
    inDir->normalize();

    *outDir = reflect(*inDir, halfDir);
    outDir->normalize();
}

template <typename ScalarT>
void HalfDifferenceCoordinateSystem::fromXyz(const Vec3&    inDir,
                                             const Vec3&    outDir,
                                             ScalarT*       halfTheta,
                                             ScalarT*       halfPhi,
                                             ScalarT*       diffTheta,
                                             ScalarT*       diffPhi)
{
    Vec3 halfDir = (inDir + outDir).normalized();
    SphericalCoordinateSystem::fromXyz(halfDir, halfTheta, halfPhi);

    Vec2 rotPhVec = Eigen::Rotation2D<Vec2::Scalar>(-*halfPhi) * Vec2(inDir[0], inDir[1]);
    Vec2 rotThVec = Eigen::Rotation2D<Vec2::Scalar>(*halfTheta) * Vec2(rotPhVec[0], inDir[2]);
    Vec3 diffDir(static_cast<Vec3::Scalar>(rotThVec[0]),
                 static_cast<Vec3::Scalar>(rotPhVec[1]),
                 static_cast<Vec3::Scalar>(rotThVec[1]));
    diffDir.normalize();

    SphericalCoordinateSystem::fromXyz(diffDir, diffTheta, diffPhi);
}

template <typename ScalarT>
void HalfDifferenceCoordinateSystem::fromXyz(const Vec3&    inDir,
                                             const Vec3&    outDir,
                                             ScalarT*       halfTheta,
                                             ScalarT*       diffTheta,
                                             ScalarT*       diffPhi)
{
    Vec3 halfDir = (inDir + outDir).normalized();
    ScalarT halfPhi; // halfPhi is 0 for isotorpic data.
    SphericalCoordinateSystem::fromXyz(halfDir, halfTheta, &halfPhi);

    Vec2 rotPhVec = Eigen::Rotation2D<Vec2::Scalar>(-halfPhi)
                  * Vec2(static_cast<Vec2::Scalar>(inDir[0]), static_cast<Vec2::Scalar>(inDir[1]));
    Vec2 rotThVec = Eigen::Rotation2D<Vec2::Scalar>(*halfTheta) * Vec2(rotPhVec[0], inDir[2]);
    Vec3 diffDir(static_cast<Vec3::Scalar>(rotThVec[0]),
                 static_cast<Vec3::Scalar>(rotPhVec[1]),
                 static_cast<Vec3::Scalar>(rotThVec[1]));
    diffDir.normalize();

    SphericalCoordinateSystem::fromXyz(diffDir, diffTheta, diffPhi);
}

} // namespace lb

#endif // LIBBSDF_HALF_DIFFERENCE_COORDINATE_SYSTEM_H
