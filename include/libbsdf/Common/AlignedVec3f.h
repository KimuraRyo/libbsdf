// =================================================================== //
// Copyright (C) 2014-2019 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#ifndef LIBBSDF_ALIGNED_VEC3F_H
#define LIBBSDF_ALIGNED_VEC3F_H

#include <Eigen/Core>
#include <Eigen/Geometry>

namespace lb {

/*!
 * \struct  AlignedVec3f
 * \brief   The AlignedVec3f struct provides an aligned vector with three components.
 */
struct AlignedVec3f : public Eigen::Vector4f
{
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    AlignedVec3f();

    AlignedVec3f(const float& x, const float& y, const float& z);

    AlignedVec3f(const Eigen::Vector3f& v);

    template<typename OtherDerived>
    AlignedVec3f(const Eigen::MatrixBase<OtherDerived>& other);

    template<typename OtherDerived>
    AlignedVec3f& operator=(const Eigen::MatrixBase<OtherDerived>& other);

    AlignedVec3f& operator=(const Eigen::Vector3f& other);

    Index rows() const;

    /*! Returns Eigen::Vector3f. */
    Eigen::Vector3f asVector3f() const;

    /*! Returns the cross product. */
    AlignedVec3f cross(const AlignedVec3f& other) const;
};

inline AlignedVec3f::AlignedVec3f() : Eigen::Vector4f() {}

inline AlignedVec3f::AlignedVec3f(const float& x, const float& y, const float& z) : Eigen::Vector4f(x, y, z, 0.0f) {}

inline AlignedVec3f::AlignedVec3f(const Eigen::Vector3f& v) : Eigen::Vector4f(v[0], v[1], v[2], 0.0f) {}

template<typename OtherDerived>
inline AlignedVec3f::AlignedVec3f(const Eigen::MatrixBase<OtherDerived>& other) : Eigen::Vector4f(other)
{
    this->data()[3] = 0.0f;
}

template<typename OtherDerived>
inline AlignedVec3f& AlignedVec3f::operator=(const Eigen::MatrixBase<OtherDerived>& other)
{
    this->Eigen::Vector4f::operator=(other);
    this->data()[3] = 0.0f;
    return *this;
}

inline AlignedVec3f& AlignedVec3f::operator=(const Eigen::Vector3f& other)
{
    this->data()[0] = other.data()[0];
    this->data()[1] = other.data()[1];
    this->data()[2] = other.data()[2];
    this->data()[3] = 0.0f;
    return *this;
}

inline AlignedVec3f::Index AlignedVec3f::rows() const { return 3; }

inline Eigen::Vector3f AlignedVec3f::asVector3f() const
{
    return Eigen::Vector3f(this->data()[0], this->data()[1], this->data()[2]);
}

inline AlignedVec3f AlignedVec3f::cross(const AlignedVec3f& other) const
{
    return this->cross3(other);
}

} // namespace lb

#endif // LIBBSDF_ALIGNED_VEC3F_H
