// =================================================================== //
// Copyright (C) 2014-2020 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#ifndef LIBBSDF_HALF_DIFFERENCE_COORDINATE_SYSTEM_H
#define LIBBSDF_HALF_DIFFERENCE_COORDINATE_SYSTEM_H

#include <libbsdf/Common/SphericalCoordinateSystem.h>

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

    static constexpr char ANGLE0_NAME[] = "half polar angle";           /*!< This attribute holds the name of halfTheta. */
    static constexpr char ANGLE1_NAME[] = "half azimuthal angle";       /*!< This attribute holds the name of halfPhi. */
    static constexpr char ANGLE2_NAME[] = "difference polar angle";     /*!< This attribute holds the name of diffTheta. */
    static constexpr char ANGLE3_NAME[] = "difference azimuthal angle"; /*!< This attribute holds the name of diffPhi. */

    static constexpr float MIN_ANGLE0 = 0.0f; /*!< This attribute holds the minimum value of halfTheta. */
    static constexpr float MIN_ANGLE1 = 0.0f; /*!< This attribute holds the minimum value of halfPhi. */
    static constexpr float MIN_ANGLE2 = 0.0f; /*!< This attribute holds the minimum value of diffTheta. */
    static constexpr float MIN_ANGLE3 = 0.0f; /*!< This attribute holds the minimum value of diffPhi. */

    static constexpr float MAX_ANGLE0 = PI_2_F; /*!< This attribute holds the maximum value of halfTheta. */
    static constexpr float MAX_ANGLE1 = TAU_F;  /*!< This attribute holds the maximum value of halfPhi. */
    static constexpr float MAX_ANGLE2 = PI_2_F; /*!< This attribute holds the maximum value of diffTheta. */
    static constexpr float MAX_ANGLE3 = TAU_F;  /*!< This attribute holds the maximum value of diffPhi. */
};

template <typename ScalarT>
void HalfDifferenceCoordinateSystem::toXyz(ScalarT  halfTheta,
                                           ScalarT  halfPhi,
                                           ScalarT  diffTheta,
                                           ScalarT  diffPhi,
                                           Vec3*    inDir,
                                           Vec3*    outDir)
{
    // Avoid unstable results.
    diffTheta = std::min(diffTheta, ScalarT(decrease(MAX_ANGLE2)));

    Vec3 halfDir = SphericalCoordinateSystem::toXyz(halfTheta, halfPhi);
    Vec3 diffDir = SphericalCoordinateSystem::toXyz(diffTheta, diffPhi);

    Vec2 rotThVec = Eigen::Rotation2D<Vec2::Scalar>(-halfTheta) * Vec2(diffDir[0], diffDir[2]);
    Vec2 rotPhVec = Eigen::Rotation2D<Vec2::Scalar>(halfPhi) * Vec2(rotThVec[0], diffDir[1]);
    *inDir = Vec3(static_cast<Vec3::Scalar>(rotPhVec[0]),
                  static_cast<Vec3::Scalar>(rotPhVec[1]),
                  static_cast<Vec3::Scalar>(rotThVec[1]));
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
    ScalarT halfPhi; // halfPhi is 0 for isotropic data.
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
