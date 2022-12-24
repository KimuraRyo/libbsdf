// =================================================================== //
// Copyright (C) 2021-2022 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#ifndef LIBBSDF_BRDF_FITTER_H
#define LIBBSDF_BRDF_FITTER_H

#include <libbsdf/Brdf/Brdf.h>
#include <libbsdf/Common/Utility.h>
#include <libbsdf/ReflectanceModel/ReflectanceModel.h>

namespace ceres {
class Problem;
}

namespace lb {

/*!
 * \class   BrdfFitter
 * \brief   The BrdfFitter class is the base class of fitters using reflectance models.
 * 
 * Parameters are acquired using Ceres Solver.
 */
class BrdfFitter
{
public:
    /*!
     * \struct  Sample
     * \brief   The Sample struct provides the sample point of BRDF.
     */
    struct Sample
    {
        Vec3 inDir;  /*!< Incoming direction. */
        Vec3 outDir; /*!< Outgoing direction. */
        Vec3 value;  /*!< BRDF value. */
    };

    /*!
     * \class   Data
     * \brief   The Data class provides sample points for fitting.
     */
    class Data
    {
    public:
        /*!
         * Constructs samples of BRDF for fitting.
         *
         * \param numSampling   The number of samples for fitting. If 0, samples in tabular data of \a brdf are used.
         * \param maxTheta      Maximum incoming and outgoing polar angle of sample points for fitting.
         */
        explicit Data(const Brdf&         brdf,
                      int                 numSampling = 100000,
                      const Vec3::Scalar& maxTheta = Vec3::Scalar(PI_2_D));

        const std::vector<Sample>& getSamples() const { return samples_; }

    private:
        std::vector<Sample> samples_;
    };

    /*! Sets the lower and upper bounds for a parameter. */
    static bool setParameterBounds(ceres::Problem*                    problem,
                                   double*                            parameter,
                                   const ReflectanceModel::Parameter& bounds);

    /*! Computes root-mean-square error between the reflectance model and data. */
    //static Vec3::Scalar computeRmse(const ReflectanceModel& model, const Data& data);

    /*! Computes error between the reflectance model and data. */
    static Vec3::Scalar computeError(const ReflectanceModel& model, const Data& data);

    /*! Computes a logarithmically scaled value. */
    template <typename T>
    static T toLogScale(const T& value);

private:
    /*! Computes a logarithmically scaled value. */
    template <typename ValueT, typename BaseT>
    static ValueT toLogScale(const ValueT& value, const BaseT& base);
};

template <typename T>
T BrdfFitter::toLogScale(const T& value)
{
    using Scalar = typename T::Scalar;
    return toLogScale(value, Scalar(10));
}

template <typename ValueT, typename BaseT>
ValueT BrdfFitter::toLogScale(const ValueT& value, const BaseT& base)
{
    return ValueT(lb::toLogScale(value[0], base),
                  lb::toLogScale(value[1], base),
                  lb::toLogScale(value[2], base));
}

} // namespace lb

#endif // LIBBSDF_BRDF_FITTER_H
