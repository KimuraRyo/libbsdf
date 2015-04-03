// =================================================================== //
// Copyright (C) 2014-2015 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#ifndef LIBBSDF_INTEGRATOR_H
#define LIBBSDF_INTEGRATOR_H

#include <cmath>

#include <libbsdf/Brdf/Brdf.h>
#include <libbsdf/Brdf/Sampler.h>

namespace lb {

/*!
 * \class   Integrator
 * \brief   The Integrator class provides functions to calculate the reflectance of a BRDF.
 *
 * Monte Carlo integration is used.
 */
class Integrator
{
public:
    /*!
     * Constructs the integrator for BRDF.
     *
     * \param numSampling                   The number of samples of Monte Carlo integration.
     * \param usePoissonDiskDistribution    If this parameter is true, precomputed Poisson disk distribution is used.
     */
    explicit Integrator(int numSampling = 100000, bool usePoissonDiskDistribution = false);

    /*! Computes the reflectance of a BRDF at an incoming direction using precomputed outgoing directions. */
    Spectrum computeReflectance(const Brdf& brdf, const Vec3& inDir);

    /*! Computes the reflectance of a BRDF at an incoming direction. */
    static Spectrum computeReflectance(const Brdf& brdf, const Vec3& inDir, int numSampling);

private:
    /*! Initializes outgoing directions for integration. */
    void initializeOutDirs(bool usePoissonDiskDistribution = false);

    int numSampling_; /*!< The number of samples of Monte Carlo integration. */

    Eigen::Array3Xf outDirs_; /*!< The array of outgoing directions. */
};

} // namespace lb

#endif // LIBBSDF_INTEGRATOR_H
