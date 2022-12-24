// =================================================================== //
// Copyright (C) 2021-2022 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#ifndef LIBBSDF_ANISOTROPIC_GGX_FITTER_H
#define LIBBSDF_ANISOTROPIC_GGX_FITTER_H

#include <libbsdf/Brdf/Brdf.h>
#include <libbsdf/Fitter/BrdfFitter.h>
#include <libbsdf/ReflectanceModel/AnisotropicGgx.h>

namespace lb {

/*!
 * \class   AnisotropicGgxFitter
 * \brief   The AnisotropicGgxFitter class provides functions to fit the parameters of the anisotropic GGX reflectance model to BRDF.
 */
class AnisotropicGgxFitter : public BrdfFitter
{
public:
    /*!
     * Estimates the fitted parameters from \a brdf.
     *
     * \param brdf          BRDF data to be fitted.
     * \param numSampling   The number of samples for fitting. If 0, samples in tabular data of \a brdf are used.
     * \param maxTheta      Maximum incoming and outgoing polar angle of sample points for fitting.
     * \return              Reflectance model with fitted parameters.
     */
    static AnisotropicGgx estimateParameters(const Brdf&         brdf,
                                             int                 numSampling = 100000,
                                             const Vec3::Scalar& maxTheta = Vec3::Scalar(PI_2_D));

    /*!
     * Estimates the fitted parameters from \a brdf.
     *
     * \param model         Reflectance model with parameters for fitting.
     * \param brdf          BRDF data to be fitted.
     * \param numSampling   The number of samples for fitting. If 0, samples in tabular data of \a brdf are used.
     * \param maxTheta      Maximum incoming and outgoing polar angle of sample points for fitting.
     */
    static void estimateParameters(AnisotropicGgx*     model,
                                   const Brdf&         brdf,
                                   int                 numSampling = 100000,
                                   const Vec3::Scalar& maxTheta = Vec3::Scalar(PI_2_D));
};

} // namespace lb

#endif // LIBBSDF_ANISOTROPIC_GGX_FITTER_H
