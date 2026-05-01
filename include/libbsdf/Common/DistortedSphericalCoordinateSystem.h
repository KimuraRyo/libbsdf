// =================================================================== //
// Copyright (C) 2026 Kimura Ryo                                       //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#ifndef LIBBSDF_DISTORTED_SPHERICAL_COORDINATE_SYSTEM_H
#define LIBBSDF_DISTORTED_SPHERICAL_COORDINATE_SYSTEM_H

#include <libbsdf/Common/Global.h>
#include <libbsdf/Common/SphericalCoordinateSystem.h>
#include <libbsdf/Common/Vector.h>

namespace lb {

/*!
 * \struct  DistortedSphericalCoordinateSystem
 * \brief   The DistortedSphericalCoordinateSystem struct provides the functions of a distorted spherical coordinate system.
 *
 * The outgoing direction is defined in a distorted spherical coordinate system with the specular direction as the zenith.
 * The coordinate system has four angle parameters.
 *   - \a inTheta: the polar angle of an incoming direction
 *   - \a inPhi: the azimuthal angle of an incoming direction
 *   - \a distTheta: the polar angle of an outgoing direction in a distorted spherical coordinate system
 *   - \a distPhi: the azimuthal angle of an outgoing direction in a distorted spherical coordinate system
 *
 * \a dist is an abbreviation for distorted. \a inPhi is not used for isotropic BRDFs.
 */
struct DistortedSphericalCoordinateSystem
{
    /*!
     * Converts from four angles to incoming and outgoing directions and
     * assigns them to \a inDir and \a outDir.
     */
    template <typename ScalarT>
    static void toXyz(ScalarT inTheta,
                      ScalarT inPhi,
                      ScalarT distTheta,
                      ScalarT distPhi,
                      Vec3*   inDir,
                      Vec3*   outDir);

    /*!
     * Converts from incoming and outgoing directions to four angles and
     * assigns them to \a inTheta, \a inPhi, \a distTheta, and \a distPhi.
     */
    template <typename ScalarT>
    static void fromXyz(const Vec3& inDir,
                        const Vec3& outDir,
                        ScalarT*    inTheta,
                        ScalarT*    inPhi,
                        ScalarT*    distTheta,
                        ScalarT*    distPhi);

    /*!
     * Converts from incoming and outgoing directions to three angles for isotropic data and
     * assigns them to \a inTheta, \a distTheta, and \a distPhi.
     */
    template <typename ScalarT>
    static void fromXyz(const Vec3& inDir,
                        const Vec3& outDir,
                        ScalarT*    inTheta,
                        ScalarT*    distTheta,
                        ScalarT*    distPhi);

    static constexpr char ANGLE0_NAME[] = "incoming polar angle";     /*!< This attribute holds the name of inTheta. */
    static constexpr char ANGLE1_NAME[] = "incoming azimuthal angle"; /*!< This attribute holds the name of inPhi. */
    static constexpr char ANGLE2_NAME[] = "distorted polar angle";     /*!< This attribute holds the name of distTheta. */
    static constexpr char ANGLE3_NAME[] = "distorted azimuthal angle"; /*!< This attribute holds the name of distPhi. */

    static constexpr double MIN_ANGLE0 = 0; /*!< This attribute holds the minimum value of inTheta. */
    static constexpr double MIN_ANGLE1 = 0; /*!< This attribute holds the minimum value of inPhi. */
    static constexpr double MIN_ANGLE2 = 0; /*!< This attribute holds the minimum value of distTheta. */
    static constexpr double MIN_ANGLE3 = 0; /*!< This attribute holds the minimum value of distPhi. */

    static constexpr double MAX_ANGLE0 = PI_2_D; /*!< This attribute holds the maximum value of inTheta. */
    static constexpr double MAX_ANGLE1 = TAU_D;  /*!< This attribute holds the maximum value of inPhi. */
    static constexpr double MAX_ANGLE2 = PI_2_D; /*!< This attribute holds the maximum value of distTheta. */
    static constexpr double MAX_ANGLE3 = TAU_D;  /*!< This attribute holds the maximum value of distPhi. */

private:
    /*! This attribute holds a small value used for numerical stability. */
    static constexpr double EPS = 1e-10;

    /*!
     * \brief Converts from distorted spherical coordinates to Cartesian coordinates.
     * \param theta     Polar angle in a distorted spherical coordinate system.
     * \param phi       Azimuthal angle in a distorted spherical coordinate system.
     * \param zenithDir Zenith direction (the direction vector pointing toward the zenith when θ = 0).
     * \return Direction vector in Cartesian coordinates.
     */
    static Vec3 distortedSphericalToXyz(double theta, double phi, const Vec3& zenithDir)
    {
        // Position on a disk in a distorted coordinate system
        double distR = 2.0 * theta / PI_D;
        Vec2   distU(distR * std::cos(phi), distR * std::sin(phi));

        // Position on the disk in the zenith direction
        Vec2 d0 = cartesianToDisk(zenithDir);

        // Mapping to a position on a standard disk
        Vec2 stdU = distU + d0 * (1.0 - distR);

        return diskToCartesian(stdU);
    }

    /*!
     * \brief Converts from Cartesian coordinates to distorted spherical coordinates.
     * \param dir       Direction vector in Cartesian coordinates.
     * \param zenithDir Zenith direction (the direction vector pointing toward the zenith when θ = 0).
     * \param theta     Polar angle in a distorted spherical coordinate system.
     * \param phi       Azimuthal angle in a distorted spherical coordinate system.
     */
    static void
    xyzToDistortedSpherical(const Vec3& dir, const Vec3& zenithDir, double* theta, double* phi)
    {
        Vec2 stdU = cartesianToDisk(dir);
        Vec2 d0 = cartesianToDisk(zenithDir);

        // Solve the quadratic equation Ac^2 + Bc + C = 0 to find c(distR).
        Vec2   diff = stdU - d0;
        double sqD0 = d0.squaredNorm();

        double A = 1.0 - sqD0;
        double B = -2.0 * diff.dot(d0);
        double C = -diff.squaredNorm();

        double c = 0.0;
        //if (std::abs(A) < EPS) {
        //    // Fallback (linear equation) when the zenith is extremely close to the equator.
        //    if (std::abs(B) > EPS) {
        //        c = -C / B;
        //    }
        //}
        //else {
            double discriminant = B * B - 4.0 * A * C;
            c = (-B + std::sqrt(discriminant)) / (2.0 * A);
        //}

        // Calculate the position distU on the disk in a distorted coordinate system.
        Vec2 distU = diff + c * d0;

        // Convert from position on a disk to spherical coordinates.
        *theta = c * PI_D / 2.0;
        *phi = std::atan2(distU.y(), distU.x());

        if (*phi < 0) {
            *phi += TAU_D;
        }
    }

    /*! Converts from a Cartesian coordinate direction vector to coordinates on a 2D disk. */
    static Vec2 cartesianToDisk(const Vec3& v)
    {
        double theta = std::acos(v.z());
        double phi = std::atan2(v.y(), v.x());
        double r = 2.0 * theta / PI_D;
        return Vec2(r * std::cos(phi), r * std::sin(phi));
    }

    /*! Converts from coordinates on a 2D disk to a Cartesian coordinate direction vector. */
    static Vec3 diskToCartesian(const Vec2& u)
    {
        double r = u.norm();
        double theta = r * PI_D / 2.0;
        double phi = std::atan2(u.y(), u.x());
        return SphericalCoordinateSystem::toXyz(theta, phi);
    }
};

template <typename ScalarT>
void DistortedSphericalCoordinateSystem::toXyz(ScalarT inTheta,
                                               ScalarT inPhi,
                                               ScalarT distTheta,
                                               ScalarT distPhi,
                                               Vec3*   inDir,
                                               Vec3*   outDir)
{
    *inDir = SphericalCoordinateSystem::toXyz(inTheta, inPhi);
    *inDir = Vec3(inDir->x(), inDir->y(), std::max(inDir->z(), EPS)).normalized();

    Vec3 reflectedInDir = Vec3(-inDir->x(), -inDir->y(), std::max(inDir->z(), EPS)).normalized();
    *outDir = distortedSphericalToXyz(distTheta, distPhi, reflectedInDir);
}

template <typename ScalarT>
void DistortedSphericalCoordinateSystem::fromXyz(const Vec3& inDir,
                                                 const Vec3& outDir,
                                                 ScalarT*    inTheta,
                                                 ScalarT*    inPhi,
                                                 ScalarT*    distTheta,
                                                 ScalarT*    distPhi)
{
    Vec3 adjustedInDir = Vec3(inDir.x(), inDir.y(), std::max(inDir.z(), EPS)).normalized();

    SphericalCoordinateSystem::fromXyz(adjustedInDir, inTheta, inPhi);

    Vec3 reflectedInDir =
        Vec3(-adjustedInDir.x(), -adjustedInDir.y(), std::max(adjustedInDir.z(), EPS)).normalized();
    xyzToDistortedSpherical(outDir, reflectedInDir, distTheta, distPhi);
}

template <typename ScalarT>
void DistortedSphericalCoordinateSystem::fromXyz(const Vec3& inDir,
                                                 const Vec3& outDir,
                                                 ScalarT*    inTheta,
                                                 ScalarT*    distTheta,
                                                 ScalarT*    distPhi)
{
    Vec3 adjustedInDir = Vec3(inDir.x(), inDir.y(), std::max(inDir.z(), EPS)).normalized();

    ScalarT inPhi;
    SphericalCoordinateSystem::fromXyz(adjustedInDir, inTheta, &inPhi);

    Vec3 reflectedInDir =
        Vec3(-adjustedInDir.x(), -adjustedInDir.y(), std::max(adjustedInDir.z(), EPS)).normalized();
    xyzToDistortedSpherical(outDir, reflectedInDir, distTheta, distPhi);

    // Convert distPhi to an angle relative to the incident direction.
    *distPhi = *distPhi - inPhi;
    if (*distPhi < 0) {
        *distPhi += static_cast<ScalarT>(TAU_D);
    }
}

} // namespace lb

#endif // LIBBSDF_DISTORTED_SPHERICAL_COORDINATE_SYSTEM_H