// =================================================================== //
// Copyright (C) 2014-2019 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

/*!
 * \file    Array.h
 * \brief   The Array.h header file includes the array declarations and functions.
 */

#ifndef LIBBSDF_ARRAY_H
#define LIBBSDF_ARRAY_H

#include <iterator>
#include <vector>

#include <Eigen/Core>

#include <libbsdf/Common/CentripetalCatmullRomSpline.h>
#include <libbsdf/Common/Utility.h>

namespace lb {

typedef Eigen::ArrayXf Arrayf;
typedef Eigen::ArrayXd Arrayd;

/*! \brief Copies an array. */
template <typename SrcT, typename DestT>
void copyArray(const SrcT& srcArray, DestT* destArray);

/*! \brief Appends an element to the end of an array. */
template <typename ArrayT, typename ScalarT>
void appendElement(ArrayT* arrayf, ScalarT value);

/*!
 * \brief Interpolates arrays using centripetal Catmull-Rom spline at \a pos in [\a pos1,\a pos2].
 * \param array Interpolated array.
 */
template <typename T>
void catmullRomSpline(float pos0, float pos1, float pos2, float pos3,
                      const T& array0, const T& array1, const T& array2, const T& array3,
                      float pos, T* array);

/*! \brief Converts an array from degrees to radians. */
template <typename T>
T toRadians(const T& degrees);

/*! \brief Returns true if the elements of an array are equally-spaced intervals. */
template <typename T>
bool isEqualInterval(const T& array);

/*!
 * Finds neighbor indices and values.
 *
 * \param lowerIndex Found index of the sample point at the lower bound.
 * \param upperIndex Found index of the sample point at the upper bound.
 * \param lowerValue Found value of the sample point at the lower bound.
 * \param upperValue Found value of the sample point at the upper bound.
 */
void findBounds(const Arrayf&   values,
                float           value,
                bool            equalIntervalValues,
                int*            lowerIndex,
                int*            upperIndex,
                float*          lowerValue,
                float*          upperValue);

/*
 * Implementation
 */

template <typename SrcT, typename DestT>
void copyArray(const SrcT& srcArray, DestT* destArray)
{
    int i = 0;
    for (auto it = srcArray.begin(); it != srcArray.end(); ++it, ++i) {
        (*destArray)[i] = *it;
    }
}

template <typename ArrayT, typename ScalarT>
void appendElement(ArrayT* arrayf, ScalarT value)
{
    ArrayT& a = *arrayf;
    std::vector<ScalarT> orig(a.data(), a.data() + a.size());
    orig.push_back(value);
    a.resize(a.size() + 1);

#if (_MSC_VER >= 1600) // Visual Studio 2010
    std::copy(orig.begin(), orig.end(),
              stdext::checked_array_iterator<ArrayT::Scalar*>(a.data(), a.size()));
#else
    std::copy(orig.begin(), orig.end(), a.data());
#endif
}

template <typename T>
void catmullRomSpline(float pos0, float pos1, float pos2, float pos3,
                      const T& array0, const T& array1, const T& array2, const T& array3,
                      float pos, T* array)
{
    assert(array0.size() == array1.size() &&
           array1.size() == array2.size() &&
           array2.size() == array3.size());
    
    array->resize(array0.size());

    CentripetalCatmullRomSpline ccrs;
    #pragma omp parallel for private(ccrs)
    for (int i = 0; i < array->size(); ++i) {
        ccrs.initialize(Vec2(pos0, array0[i]),
                        Vec2(pos1, array1[i]),
                        Vec2(pos2, array2[i]),
                        Vec2(pos3, array3[i]));

        typedef typename T::Scalar ScalarType;
        (*array)[i] = static_cast<ScalarType>(ccrs.interpolateY(pos));
    }
}

template <typename T>
T toRadians(const T& degrees)
{
    typedef typename T::Scalar ScalarType;
    return degrees / ScalarType(180) * ScalarType(PI_D);
}

template <typename T>
bool isEqualInterval(const T& array)
{
    if (array.size() <= 2) return false;

    float interval = array[array.size() - 1] / (array.size() - 1);
    for (int i = 0; i < array.size(); ++i) {
        if (!isEqual(array[i], interval * i)) {
            return false;
        }
    }

    return true;
}

} // namespace lb

#endif // LIBBSDF_ARRAY_H
