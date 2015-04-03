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

namespace lb {

/*!
 * \class   SpecularCoordinateSystem
 * \brief   The SpecularCoordinateSystem class provides the functions of a specular coordinate system.
 *
 * The coordinate system has four angle parameters.
 *   - inTheta: the polar angle of an incoming direction
 *   - inPhi: the azimuthal angle of an incoming direction
 *   - specTheta: the polar angle of a specular direction
 *   - specPhi: the azimuthal angle of a specular direction
 *
 * Spec is an abbreviation for specular.
 */
class SpecularCoordinateSystem
{
public:
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
    static Vec3 toOutDirXyz(float inTheta, float specTheta, float specPhi);

    /*! Converts an outgoing direction from a Cartesian coordinate system to a specular. */
    static void fromOutDirXyz(const Vec3& dir, float inTheta,
                              float* specTheta, float* specPhi);
};

inline void SpecularCoordinateSystem::toXyz(float inTheta, float inPhi,
                                            float specTheta, float specPhi,
                                            Vec3* inDir, Vec3* outDir)
{
    *inDir = SphericalCoordinateSystem::toXyz(inTheta, inPhi);
    *outDir = toOutDirXyz(inTheta, specTheta, specPhi);
}

inline void SpecularCoordinateSystem::fromXyz(const Vec3& inDir, const Vec3& outDir,
                                              float* inTheta, float* inPhi,
                                              float* specTheta, float* specPhi)
{
    SphericalCoordinateSystem::fromXyz(inDir, inTheta, inPhi);
    fromOutDirXyz(outDir, *inTheta, specTheta, specPhi);
}

inline Vec3 SpecularCoordinateSystem::toOutDirXyz(float inTheta, float specTheta, float specPhi)
{
    Vec3 xyzVec = SphericalCoordinateSystem::toXyz(specTheta, specPhi);
    Vec2f rotVec = Eigen::Rotation2D<Vec2f::Scalar>(inTheta) * Vec2f(xyzVec(0), xyzVec(2));

    return Vec3(rotVec(0), xyzVec(1), rotVec(1));
}

inline void SpecularCoordinateSystem::fromOutDirXyz(const Vec3& dir, float inTheta,
                                                    float* specTheta, float* specPhi)
{
    Vec2f rotVec = Eigen::Rotation2D<Vec2f::Scalar>(-inTheta) * Vec2f(dir(0), dir(2));
    rotVec(1) = std::min(rotVec(1), 1.0f);
    *specTheta = std::acos(rotVec(1));
    *specPhi = std::atan2(dir(1), rotVec(0));
    if (*specPhi < 0.0f) {
        *specPhi += 2.0f * PI_F;
    }
}

} // namespace lb

#endif // LIBBSDF_SPECULAR_COORDINATE_SYSTEM_H
