// =================================================================== //
// Copyright (C) 2021-2022 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#ifndef LIBBSDF_UNREAL_ENGINE_4_FITTER_H
#define LIBBSDF_UNREAL_ENGINE_4_FITTER_H

#include <libbsdf/Brdf/Brdf.h>
#include <libbsdf/Fitter/BrdfFitter.h>
#include <libbsdf/ReflectanceModel/UnrealEngine4.h>

namespace lb {

/*!
 * \class   UnrealEngine4Fitter
 * \brief   The UnrealEngine4Fitter class provides functions to fit the parameters of the Unreal Engine 4 reflectance model to BRDF.
 */
class UnrealEngine4Fitter : public BrdfFitter
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
    static UnrealEngine4 estimateParameters(const Brdf&         brdf,
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
    static void estimateParameters(UnrealEngine4*      model,
                                   const Brdf&         brdf,
                                   int                 numSampling = 100000,
                                   const Vec3::Scalar& maxTheta = Vec3::Scalar(PI_2_D));

private:
    /*! Acquires fitted parameters using a constant metallic value. This function is used internally for better fitting. */
    static void estimateWithoutMetallic(UnrealEngine4* model, const Data& data);
};

} // namespace lb

#endif // LIBBSDF_UNREAL_ENGINE_4_FITTER_H
