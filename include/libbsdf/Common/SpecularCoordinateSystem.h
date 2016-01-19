// =================================================================== //
// Copyright (C) 2014-2015 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#ifndef LIBBSDF_SPECULAR_COORDINATE_SYSTEM_H
#define LIBBSDF_SPECULAR_COORDINATE_SYSTEM_H

#include <Eigen/Geometry>

#include <libbsdf/Common/SphericalCoordinateSystem.h>
#include <libbsdf/Common/Utility.h>

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
 * \a spec is an abbreviation for specular. \a inPhi isn't used for isotropic BRDFs.
 */
struct SpecularCoordinateSystem
{
    /*!
     * Converts from four angles to incoming and outgoing directions and
     * assigns them to \a inDir and \a outDir.
     */
    static void toXyz(float inTheta, float inPhi,
                      float specTheta, float specPhi,
                      Vec3* inDir, Vec3* outDir);

    /*!
     * Converts from incoming and outgoing directions to four angles and
     * assigns them to \a inTheta, \a inPhi, \a specTheta, and \a specPhi.
     */
    static void fromXyz(const Vec3& inDir, const Vec3& outDir,
                        float* inTheta, float* inPhi,
                        float* specTheta, float* specPhi);

    /*!
     * Converts from incoming and outgoing directions to three angles for isotropic data and
     * assigns them to \a inTheta, \a specTheta, and \a specPhi.
     */
    static void fromXyz(const Vec3& inDir, const Vec3& outDir,
                        float* inTheta,
                        float* specTheta, float* specPhi);

    static const std::string ANGLE0_NAME; /*!< This attribute holds the name of inTheta. */
    static const std::string ANGLE1_NAME; /*!< This attribute holds the name of inPhi. */
    static const std::string ANGLE2_NAME; /*!< This attribute holds the name of specTheta. */
    static const std::string ANGLE3_NAME; /*!< This attribute holds the name of specPhi. */

    static const float MAX_ANGLE0; /*!< This attribute holds the maximum value of inTheta. */
    static const float MAX_ANGLE1; /*!< This attribute holds the maximum value of inPhi. */
    static const float MAX_ANGLE2; /*!< This attribute holds the maximum value of specTheta. */
    static const float MAX_ANGLE3; /*!< This attribute holds the maximum value of specPhi. */

private:
    /*! Converts an outgoing direction from a specular coordinate system to a Cartesian. */
    static Vec3 toOutDirXyz(float inTheta, float inPhi,
                            float specTheta, float specPhi);

    /*! Converts an outgoing direction from a Cartesian coordinate system to a specular. */
    static void fromOutDirXyz(const Vec3& outDir, float inTheta, float inPhi,
                              float* specTheta, float* specPhi);
};

inline void SpecularCoordinateSystem::toXyz(float inTheta, float inPhi,
                                            float specTheta, float specPhi,
                                            Vec3* inDir, Vec3* outDir)
{
    *inDir = SphericalCoordinateSystem::toXyz(inTheta, inPhi);
    *outDir = toOutDirXyz(inTheta, inPhi, specTheta, specPhi);
}

inline void SpecularCoordinateSystem::fromXyz(const Vec3& inDir, const Vec3& outDir,
                                              float* inTheta, float* inPhi,
                                              float* specTheta, float* specPhi)
{
    SphericalCoordinateSystem::fromXyz(inDir, inTheta, inPhi);
    fromOutDirXyz(outDir, *inTheta, *inPhi, specTheta, specPhi);
}

inline void SpecularCoordinateSystem::fromXyz(const Vec3& inDir, const Vec3& outDir,
                                              float* inTheta,
                                              float* specTheta, float* specPhi)
{
    float inPhi; // inPhi is 0 for isotorpic data.
    SphericalCoordinateSystem::fromXyz(inDir, inTheta, &inPhi);
    fromOutDirXyz(outDir, *inTheta, inPhi, specTheta, specPhi);
}

inline Vec3 SpecularCoordinateSystem::toOutDirXyz(float inTheta, float inPhi,
                                                  float specTheta, float specPhi)
{
    Vec3 outDir = SphericalCoordinateSystem::toXyz(specTheta, specPhi);
    Vec2f rotThVec = Eigen::Rotation2D<Vec2f::Scalar>(inTheta) * Vec2f(outDir[0], outDir[2]);
    Vec2f rotPhVec = Eigen::Rotation2D<Vec2f::Scalar>(inPhi) * Vec2f(rotThVec[0], outDir[1]);

    return Vec3(rotPhVec[0], rotPhVec[1], rotThVec[1]);
}

inline void SpecularCoordinateSystem::fromOutDirXyz(const Vec3& outDir, float inTheta, float inPhi,
                                                    float* specTheta, float* specPhi)
{
    Vec2f rotPhVec = Eigen::Rotation2D<Vec2f::Scalar>(-inPhi) * Vec2f(outDir[0], outDir[1]);
    Vec2f rotThVec = Eigen::Rotation2D<Vec2f::Scalar>(-inTheta) * Vec2f(rotPhVec[0], outDir[2]);
    rotThVec[1] = clamp(rotThVec[1], -1.0f, 1.0f);
    Vec3 rotDir(rotThVec[0], rotPhVec[1], rotThVec[1]);
    SphericalCoordinateSystem::fromXyz(rotDir, specTheta, specPhi);
}

} // namespace lb

#endif // LIBBSDF_SPECULAR_COORDINATE_SYSTEM_H
