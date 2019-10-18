// =================================================================== //
// Copyright (C) 2014-2019 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#ifndef LIBBSDF_SPHERICAL_COORDINATE_SYSTEM_H
#define LIBBSDF_SPHERICAL_COORDINATE_SYSTEM_H

#include <libbsdf/Common/Global.h>
#include <libbsdf/Common/Utility.h>
#include <libbsdf/Common/Vector.h>

namespace lb {

/*!
 * \struct  SphericalCoordinateSystem
 * \brief   The SphericalCoordinateSystem struct provides the functions of a spherical coordinate system.
 *
 * The coordinate system has four angle parameters.
 *   - \a inTheta: the polar angle of an incoming direction
 *   - \a inPhi: the azimuthal angle of an incoming direction
 *   - \a outTheta: the polar angle of an outgoing direction
 *   - \a outPhi: the azimuthal angle of an outgoing direction
 *
 * \a inPhi is not used for isotropic BRDFs.
 */
struct SphericalCoordinateSystem
{
    /*!
     * Converts from four angles to incoming and outgoing directions and
     * assigns them to \a inDir and \a outDir.
     */
    template <typename ScalarT>
    static void toXyz(ScalarT   inTheta,
                      ScalarT   inPhi,
                      ScalarT   outTheta,
                      ScalarT   outPhi,
                      Vec3*     inDir,
                      Vec3*     outDir);

    /*!
     * Converts from incoming and outgoing directions to four angles and
     * assigns them to \a inTheta, \a inPhi, \a outTheta, and \a outPhi.
     */
    template <typename ScalarT>
    static void fromXyz(const Vec3& inDir,
                        const Vec3& outDir,
                        ScalarT*    inTheta,
                        ScalarT*    inPhi,
                        ScalarT*    outTheta,
                        ScalarT*    outPhi);

    /*!
     * Converts from incoming and outgoing directions to three angles for isotropic data and
     * assigns them to \a inTheta, \a outTheta, and \a outPhi.
     */
    template <typename ScalarT>
    static void fromXyz(const Vec3& inDir,
                        const Vec3& outDir,
                        ScalarT*    inTheta,
                        ScalarT*    outTheta,
                        ScalarT*    outPhi);

    static constexpr char ANGLE0_NAME[] = "Incoming polar angle";     /*!< This attribute holds the name of inTheta. */
    static constexpr char ANGLE1_NAME[] = "Incoming azimuthal angle"; /*!< This attribute holds the name of inPhi. */
    static constexpr char ANGLE2_NAME[] = "Outgoing polar angle";     /*!< This attribute holds the name of outTheta. */
    static constexpr char ANGLE3_NAME[] = "Outgoing azimuthal angle"; /*!< This attribute holds the name of outPhi. */

    static constexpr float MIN_ANGLE0 = 0.0f; /*!< This attribute holds the minimum value of inTheta. */
    static constexpr float MIN_ANGLE1 = 0.0f; /*!< This attribute holds the minimum value of inPhi. */
    static constexpr float MIN_ANGLE2 = 0.0f; /*!< This attribute holds the minimum value of outTheta. */
    static constexpr float MIN_ANGLE3 = 0.0f; /*!< This attribute holds the minimum value of outPhi. */

    static constexpr float MAX_ANGLE0 = decrease(PI_2_F);   /*!< This attribute holds the maximum value of inTheta. */
    static constexpr float MAX_ANGLE1 = decrease(TAU_F);    /*!< This attribute holds the maximum value of inPhi. */
    static constexpr float MAX_ANGLE2 = decrease(PI_2_F);   /*!< This attribute holds the maximum value of outTheta. */
    static constexpr float MAX_ANGLE3 = decrease(TAU_F);    /*!< This attribute holds the maximum value of outPhi. */

    /*! Converts from a spherical coordinate system to a Cartesian. */
    template <typename ScalarT>
    static Vec3 toXyz(ScalarT theta, ScalarT phi);

    /*! Converts from a Cartesian coordinate system to a spherical. */
    template <typename ScalarT>
    static void fromXyz(const Vec3& dir, ScalarT* theta, ScalarT* phi);
};

template <typename ScalarT>
void SphericalCoordinateSystem::toXyz(ScalarT   inTheta,
                                      ScalarT   inPhi,
                                      ScalarT   outTheta,
                                      ScalarT   outPhi,
                                      Vec3*     inDir,
                                      Vec3*     outDir)
{
    Vec4 thetaPhi(inTheta, inPhi, outTheta, outPhi);
    Vec4 sinArray = thetaPhi.array().sin();
    Vec4 cosArray = thetaPhi.array().cos();

    Vec4 lhs(sinArray[0], sinArray[0], sinArray[2], sinArray[2]);
    Vec4 rhs(cosArray[1], sinArray[1], cosArray[3], sinArray[3]);
    Vec4 productArray = lhs.cwiseProduct(rhs);

    *inDir  = Vec3(productArray[0], productArray[1], cosArray[0]).normalized();
    *outDir = Vec3(productArray[2], productArray[3], cosArray[2]).normalized();
}

template <typename ScalarT>
void SphericalCoordinateSystem::fromXyz(const Vec3& inDir,
                                        const Vec3& outDir,
                                        ScalarT*    inTheta,
                                        ScalarT*    inPhi,
                                        ScalarT*    outTheta,
                                        ScalarT*    outPhi)
{
    fromXyz(inDir, inTheta, inPhi);
    fromXyz(outDir, outTheta, outPhi);
}

template <typename ScalarT>
void SphericalCoordinateSystem::fromXyz(const Vec3& inDir,
                                        const Vec3& outDir,
                                        ScalarT*    inTheta,
                                        ScalarT*    outTheta,
                                        ScalarT*    outPhi)
{
    ScalarT inPhi; // inPhi is 0 for isotorpic data.
    fromXyz(inDir, inTheta, &inPhi);
    fromXyz(outDir, outTheta, outPhi);

    *outPhi = *outPhi - inPhi;
    if (*outPhi < 0.0) {
        *outPhi += ScalarT(TAU_D);
    }
}

template <typename ScalarT>
Vec3 SphericalCoordinateSystem::toXyz(ScalarT theta, ScalarT phi)
{
    Vec2f thetaPhi(theta, phi);
    Vec2f sinArray = thetaPhi.array().sin();
    Vec2f cosArray = thetaPhi.array().cos();

    Vec3 lhs(sinArray[0], sinArray[0], cosArray[0]);
    Vec3 rhs(cosArray[1], sinArray[1], 1);

    return lhs.cwiseProduct(rhs).normalized();
}

template <typename ScalarT>
void SphericalCoordinateSystem::fromXyz(const Vec3& dir, ScalarT* theta, ScalarT* phi)
{
    *theta = static_cast<ScalarT>(std::acos(dir[2]));
    *phi   = static_cast<ScalarT>(std::atan2(dir[1], dir[0]));
    if (*phi < 0.0) {
        *phi += ScalarT(TAU_D);
    }

    assert(!std::isnan(*theta) && !std::isnan(*phi) &&
           !std::isinf(*theta) && !std::isinf(*phi));
}

} // namespace lb

#endif // LIBBSDF_SPHERICAL_COORDINATE_SYSTEM_H
