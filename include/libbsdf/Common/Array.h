// =================================================================== //
// Copyright (C) 2014-2015 Kimura Ryo                                  //
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

#include <vector>

#include <Eigen/Core>

#include <libbsdf/Common/Utility.h>

namespace lb {

typedef Eigen::ArrayXf Arrayf;
typedef Eigen::ArrayXd Arrayd;

/*! \brief Copies an array. */
template <typename SrcT, typename DestT>
void copyArray(const SrcT& srcArray, DestT& destArray, int size);

/*! \brief Copies an array. */
template <typename SrcT, typename DestT>
void copyArray(const SrcT& srcArray, DestT& destArray);

/*! \brief Appends an element to the end of an array. */
template <typename ArrayT, typename ScalarT>
void appendElement(ArrayT& arrayf, ScalarT value);

/*! \brief Converts an array from degree to radian. */
template <typename T>
T toRadians(const T& degrees);

/*! \brief Returns true if the elements of an array are equally-spaced intervals. */
template <typename T>
bool isEqualInterval(const T& array);

/*
 * Implementation
 */

template <typename SrcT, typename DestT>
inline void copyArray(const SrcT& srcArray, DestT& destArray, int size)
{
    for (int i = 0; i < size; ++i) {
        destArray[i] = srcArray[i];
    }
}

template <typename SrcT, typename DestT>
inline void copyArray(const SrcT& srcArray, DestT& destArray)
{
    int i = 0;
    for (auto it = srcArray.begin(); it != srcArray.end(); ++it, ++i) {
        destArray[i] = *it;
    }
}

template <typename ArrayT, typename ScalarT>
inline void appendElement(ArrayT& arrayf, ScalarT value)
{
    std::vector<ScalarT> orig(arrayf.data(), arrayf.data() + arrayf.size());
    orig.push_back(value);
    arrayf.resize(arrayf.size() + 1);

#if (_MSC_VER >= 1600) // Visual Studio 2010
    std::copy(orig.begin(), orig.end(),
              stdext::checked_array_iterator<ArrayT::Scalar*>(arrayf.data(), arrayf.size()));
#else
    std::copy(orig.begin(), orig.end(), arrayf.data());
#endif
}

template <typename T>
inline T toRadians(const T& degrees)
{
    typedef typename T::Scalar ScalarType;
    return degrees / static_cast<ScalarType>(180.0) * static_cast<ScalarType>(PI_F);
}

template <typename T>
inline bool isEqualInterval(const T& array)
{
    if (array.size() <= 1) return false;

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
