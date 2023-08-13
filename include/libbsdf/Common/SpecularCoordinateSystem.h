// =================================================================== //
// Copyright (C) 2014-2023 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#ifndef LIBBSDF_SPECULAR_COORDINATE_SYSTEM_H
#define LIBBSDF_SPECULAR_COORDINATE_SYSTEM_H

#include <libbsdf/Common/SphericalCoordinateSystem.h>

namespace lb {

/*!
 * \struct  SpecularCoordinateSystem
 * \brief   The SpecularCoordinateSystem struct provides the functions of a specular coordinate system.
 *
 * The coordinate system has four angle parameters.
 *   - \a inTheta: the polar angle of an incoming direction
 *   - \a inPhi: the azimuthal angle of an incoming direction
 *   - \a specTheta: the polar angle of a specular direction
 *   - \a specPhi: the azimuthal angle of a specular direction
 *
 * \a spec is an abbreviation for specular. \a inPhi is not used for isotropic BRDFs.
 */
struct SpecularCoordinateSystem
{
    /*!
     * Converts from four angles to incoming and outgoing directions and
     * assigns them to \a inDir and \a outDir.
     */
    template <typename ScalarT>
    static void toXyz(ScalarT inTheta,
                      ScalarT inPhi,
                      ScalarT specTheta,
                      ScalarT specPhi,
                      Vec3*   inDir,
                      Vec3*   outDir);

    /*!
     * Converts from incoming and outgoing directions to four angles and
     * assigns them to \a inTheta, \a inPhi, \a specTheta, and \a specPhi.
     */
    template <typename ScalarT>
    static void fromXyz(const Vec3& inDir,
                        const Vec3& outDir,
                        ScalarT*    inTheta,
                        ScalarT*    inPhi,
                        ScalarT*    specTheta,
                        ScalarT*    specPhi);

    /*!
     * Converts from incoming and outgoing directions to three angles for isotropic data and
     * assigns them to \a inTheta, \a specTheta, and \a specPhi.
     */
    template <typename ScalarT>
    static void fromXyz(const Vec3& inDir,
                        const Vec3& outDir,
                        ScalarT*    inTheta,
                        ScalarT*    specTheta,
                        ScalarT*    specPhi);

    static constexpr char ANGLE0_NAME[] = "incoming polar angle";     /*!< This attribute holds the name of inTheta. */
    static constexpr char ANGLE1_NAME[] = "incoming azimuthal angle"; /*!< This attribute holds the name of inPhi. */
    static constexpr char ANGLE2_NAME[] = "specular polar angle";     /*!< This attribute holds the name of specTheta. */
    static constexpr char ANGLE3_NAME[] = "specular azimuthal angle"; /*!< This attribute holds the name of specPhi. */

    static constexpr double MIN_ANGLE0 = 0; /*!< This attribute holds the minimum value of inTheta. */
    static constexpr double MIN_ANGLE1 = 0; /*!< This attribute holds the minimum value of inPhi. */
    static constexpr double MIN_ANGLE2 = 0; /*!< This attribute holds the minimum value of specTheta. */
    static constexpr double MIN_ANGLE3 = 0; /*!< This attribute holds the minimum value of specPhi. */

    static constexpr double MAX_ANGLE0 = PI_2_D; /*!< This attribute holds the maximum value of inTheta. */
    static constexpr double MAX_ANGLE1 = TAU_D;  /*!< This attribute holds the maximum value of inPhi. */
    static constexpr double MAX_ANGLE2 = PI_D;   /*!< This attribute holds the maximum value of specTheta. */
    static constexpr double MAX_ANGLE3 = TAU_D;  /*!< This attribute holds the maximum value of specPhi. */

    /*! Converts an outgoing direction from a specular coordinate system to a Cartesian. */
    template <typename ScalarT>
    static Vec3 toOutDirXyz(ScalarT inTheta, ScalarT inPhi, ScalarT specTheta, ScalarT specPhi);

    /*! Converts an outgoing direction from a Cartesian coordinate system to a specular. */
    template <typename ScalarT>
    static void fromOutDirXyz(const Vec3& outDir,
                              ScalarT     inTheta,
                              ScalarT     inPhi,
                              ScalarT*    specTheta,
                              ScalarT*    specPhi);
};

template <typename ScalarT>
void SpecularCoordinateSystem::toXyz(ScalarT inTheta,
                                     ScalarT inPhi,
                                     ScalarT specTheta,
                                     ScalarT specPhi,
                                     Vec3*   inDir,
                                     Vec3*   outDir)
{
    *inDir = SphericalCoordinateSystem::toXyz(inTheta, inPhi);
    *outDir = toOutDirXyz(inTheta, inPhi, specTheta, specPhi);
}

template <typename ScalarT>
void SpecularCoordinateSystem::fromXyz(const Vec3& inDir,
                                       const Vec3& outDir,
                                       ScalarT*    inTheta,
                                       ScalarT*    inPhi,
                                       ScalarT*    specTheta,
                                       ScalarT*    specPhi)
{
    SphericalCoordinateSystem::fromXyz(inDir, inTheta, inPhi);
    fromOutDirXyz(outDir, *inTheta, *inPhi, specTheta, specPhi);
}

template <typename ScalarT>
void SpecularCoordinateSystem::fromXyz(const Vec3& inDir,
                                       const Vec3& outDir,
                                       ScalarT*    inTheta,
                                       ScalarT*    specTheta,
                                       ScalarT*    specPhi)
{
    ScalarT inPhi; // inPhi is 0 for isotropic data.
    SphericalCoordinateSystem::fromXyz(inDir, inTheta, &inPhi);
    fromOutDirXyz(outDir, *inTheta, inPhi, specTheta, specPhi);
}

template <typename ScalarT>
Vec3 SpecularCoordinateSystem::toOutDirXyz(ScalarT inTheta,
                                           ScalarT inPhi,
                                           ScalarT specTheta,
                                           ScalarT specPhi)
{
    Vec3 outDir = SphericalCoordinateSystem::toXyz(specTheta, specPhi);
    Vec2 rotThVec = Eigen::Rotation2D<Vec2::Scalar>(inTheta) * Vec2(outDir[0], outDir[2]);
    Vec2 rotPhVec = Eigen::Rotation2D<Vec2::Scalar>(inPhi) * Vec2(rotThVec[0], outDir[1]);

    Vec3 rotDir = Vec3(static_cast<Vec3::Scalar>(rotPhVec[0]),
                       static_cast<Vec3::Scalar>(rotPhVec[1]),
                       static_cast<Vec3::Scalar>(rotThVec[1]));

    return rotDir.normalized();
}

template <typename ScalarT>
void SpecularCoordinateSystem::fromOutDirXyz(const Vec3& outDir,
                                             ScalarT     inTheta,
                                             ScalarT     inPhi,
                                             ScalarT*    specTheta,
                                             ScalarT*    specPhi)
{
    Vec2 rotPhVec = Eigen::Rotation2D<Vec2::Scalar>(-inPhi) * Vec2(outDir[0], outDir[1]);
    Vec2 rotThVec = Eigen::Rotation2D<Vec2::Scalar>(-inTheta) * Vec2(rotPhVec[0], outDir[2]);
    rotThVec[1] = clamp(rotThVec[1], Vec2::Scalar(-1), Vec2::Scalar(1));

    Vec3 rotDir(static_cast<Vec3::Scalar>(rotThVec[0]),
                static_cast<Vec3::Scalar>(rotPhVec[1]),
                static_cast<Vec3::Scalar>(rotThVec[1]));
    rotDir.normalize();

    SphericalCoordinateSystem::fromXyz(rotDir, specTheta, specPhi);
}

} // namespace lb

#endif // LIBBSDF_SPECULAR_COORDINATE_SYSTEM_H
