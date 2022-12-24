// =================================================================== //
// Copyright (C) 2022 Kimura Ryo                                       //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#ifndef LIBBSDF_SIMPLE_GGX_FITTER_H
#define LIBBSDF_SIMPLE_GGX_FITTER_H

#include <libbsdf/Brdf/Brdf.h>
#include <libbsdf/Fitter/BrdfFitter.h>
#include <libbsdf/ReflectanceModel/SimpleGgx.h>

namespace lb {

/*!
 * \class   SimpleGgxFitter
 * \brief   The SimpleGgxFitter class provides functions to fit the parameters of the GGX reflectance model to BRDF.
 */
class SimpleGgxFitter : public BrdfFitter
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
    static SimpleGgx estimateParameters(const Brdf&         brdf,
                                        int                 numSampling = 100000,
                                        const Vec3::Scalar& maxTheta = Vec3::Scalar(PI_2_D));

    /*!
     * Estimates the fitted parameters from \a brdf.
     *
     * \attention   Initial values of parameters of \a model can affect fitting results.
     *
     * \param model         Reflectance model with parameters for fitting.
     * \param brdf          BRDF data to be fitted.
     * \param numSampling   The number of samples for fitting. If 0, samples in tabular data of \a brdf are used.
     * \param maxTheta      Maximum incoming and outgoing polar angle of sample points for fitting.
     */
    static void estimateParameters(SimpleGgx*          model,
                                   const Brdf&         brdf,
                                   int                 numSampling = 100000,
                                   const Vec3::Scalar& maxTheta = Vec3::Scalar(PI_2_D));
};

} // namespace lb

#endif // LIBBSDF_SIMPLE_GGX_FITTER_H
