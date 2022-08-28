// =================================================================== //
// Copyright (C) 2021 Kimura Ryo                                       //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#ifndef LIBBSDF_LAMBERTIAN_FITTER_H
#define LIBBSDF_LAMBERTIAN_FITTER_H

#include <libbsdf/Brdf/Brdf.h>
#include <libbsdf/Fitter/BrdfFitter.h>
#include <libbsdf/ReflectanceModel/Lambertian.h>

namespace lb {

/*!
 * \class   LambertianFitter
 * \brief   The LambertianFitter class provides functions to fit the parameters of the Lambertian reflectance model to BRDF.
 */
class LambertianFitter : public BrdfFitter
{
public:
    /*!
     * Estimates the fitted parameters from \a brdf.
     * 
     * \param model         Reflectance model including the parameters.
     * \param brdf          BRDF data to be fitted.
     * \param numSampling   The number of samples for fitting. If 0, samples in tabular data of lb::Brdf is used.
     * \param maxTheta      Maximum polar angle to avoid using data from large polar angles.
     */
    static void estimateParameters(Lambertian*         model,
                                   const Brdf&         brdf,
                                   int                 numSampling = 100000,
                                   const Vec3::Scalar& maxTheta = Vec3::Scalar(PI_2_D));
};

} // namespace lb

#endif // LIBBSDF_LAMBERTIAN_FITTER_H
