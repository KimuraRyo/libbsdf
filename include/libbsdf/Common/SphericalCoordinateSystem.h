// =================================================================== //
// Copyright (C) 2014-2015 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#ifndef LIBBSDF_SPHERICAL_COORDINATE_SYSTEM_H
#define LIBBSDF_SPHERICAL_COORDINATE_SYSTEM_H

#include <libbsdf/Common/Global.h>
#include <libbsdf/Common/Vector.h>

namespace lb {

/*!
 * \class   SphericalCoordinateSystem
 * \brief   The SphericalCoordinateSystem class provides the functions of a spherical coordinate system.
 *
 * The coordinate system has four angle parameters.
 *   - inTheta: the polar angle of an incoming direction
 *   - inPhi: the azimuthal angle of an incoming direction
 *   - outTheta: the polar angle of an outgoing direction
 *   - outPhi: the azimuthal angle of an outgoing direction
 */
class SphericalCoordinateSystem
{
public:
    /*!
     * Converts from four angles to incoming and outgoing directions and
     * assigns them to *inDir and *outDir.
     */
    static void toXyz(float inTheta, float inPhi,
                      float outTheta, float outPhi,
                      Vec3* inDir, Vec3* outDir);

    /*!
     * Converts from incoming and outgoing directions to four angles and
     * assigns them to *inTheta, *inPhi, *outTheta, and *outPhi.
     */
    static void fromXyz(const Vec3& inDir, const Vec3& outDir,
                        float* inTheta, float* inPhi,
                        float* outTheta, float* outPhi);

    static const std::string ANGLE0_NAME; /*!< This attribute holds the name of inTheta. */
    static const std::string ANGLE1_NAME; /*!< This attribute holds the name of inPhi. */
    static const std::string ANGLE2_NAME; /*!< This attribute holds the name of outTheta. */
    static const std::string ANGLE3_NAME; /*!< This attribute holds the name of outPhi. */

    static const float MAX_ANGLE0; /*!< This attribute holds the maximum value of inTheta. */
    static const float MAX_ANGLE1; /*!< This attribute holds the maximum value of inPhi. */
    static const float MAX_ANGLE2; /*!< This attribute holds the maximum value of outTheta. */
    static const float MAX_ANGLE3; /*!< This attribute holds the maximum value of outPhi. */

    /*! Converts from a spherical coordinate system to a Cartesian. */
    static Vec3 toXyz(float theta, float phi);

    /*! Converts from a Cartesian coordinate system to a spherical. */
    static void fromXyz(const Vec3& dir, float* theta, float* phi);
};

inline void SphericalCoordinateSystem::toXyz(float inTheta, float inPhi,
                                             float outTheta, float outPhi,
                                             Vec3* inDir, Vec3* outDir)
{
    Vec4 thetaPhi(inTheta, inPhi, outTheta, outPhi);
    Vec4 sinArray = thetaPhi.array().sin();
    Vec4 cosArray = thetaPhi.array().cos();

    Vec4 lhs(sinArray[0], sinArray[0], sinArray[2], sinArray[2]);
    Vec4 rhs(cosArray[1], sinArray[1], cosArray[3], sinArray[3]);
    Vec4 productArray = lhs.cwiseProduct(rhs);

    *inDir = Vec3(productArray[0], productArray[1], cosArray[0]);
    *outDir = Vec3(productArray[2], productArray[3], cosArray[2]);
}

inline void SphericalCoordinateSystem::fromXyz(const Vec3& inDir, const Vec3& outDir,
                                               float* inTheta, float* inPhi,
                                               float* outTheta, float* outPhi)
{
    fromXyz(inDir, inTheta, inPhi);
    fromXyz(outDir, outTheta, outPhi);
}

inline Vec3 SphericalCoordinateSystem::toXyz(float theta, float phi)
{
    Vec2f thetaPhi(theta, phi);
    Vec2f sinArray = thetaPhi.array().sin();
    Vec2f cosArray = thetaPhi.array().cos();

    Vec3 lhs(sinArray[0], sinArray[0], cosArray[0]);
    Vec3 rhs(cosArray[1], sinArray[1], 1.0f);

    return lhs.cwiseProduct(rhs);
}

inline void SphericalCoordinateSystem::fromXyz(const Vec3& dir, float* theta, float* phi)
{
    *theta = std::acos(dir[2]);
    *phi = std::atan2(dir[1], dir[0]);
    if (*phi < 0.0f) {
        *phi += 2.0f * PI_F;
    }

    assert(!std::isnan(*theta) && !std::isnan(*phi));
}

} // namespace lb

#endif // LIBBSDF_SPHERICAL_COORDINATE_SYSTEM_H
