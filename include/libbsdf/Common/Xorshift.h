// =================================================================== //
// Copyright (C) 2014-2022 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#ifndef LIBBSDF_XORSHIFT_H
#define LIBBSDF_XORSHIFT_H

#include <limits>

namespace lb {

#if defined(__C99__) || (defined(__GNUC__) && __GNUC__ >= 3)
    #include <inttypes.h>
    #include <stdint.h>
#else
    #if defined(__GNUC__)
        using uint32_t = unsigned int;
    #elif defined(_MSC_VER) || defined(__BORLANDC__)
        using uint32_t = unsigned __int32;
    #endif
#endif

/*!
 * \class   Xorshift
 * \brief   The Xorshift class provides a random number generator using Xorshift.
 */
class Xorshift
{
public:
    explicit Xorshift(uint32_t seed = 123456789);

    void setSeed(uint32_t seed);

    /*! Generates a random integer. Range is [0,std::numeric_limits<uint32_t>::max()]. */
    uint32_t next();

    /*! Generates a random floating-point number. Range is [0.0,1.0]. */
    template <typename T>
    T next();

    /*! Generates a random point on the surface of a unit hemisphere. The coordinate system is Z-up. */
    template <typename Vec3T>
    Vec3T nextOnHemisphere();

    /*! Generates a random integer. Range is [0,std::numeric_limits<uint32_t>::max()]. */
    static uint32_t random();

    /*! Generates a random floating-point number. Range is [0.0,1.0]. */
    template <typename T>
    static T random();

    /*! Generates a random point on the surface of a unit hemisphere. The coordinate system is Z-up. */
    template <typename Vec3T>
    static Vec3T randomOnHemisphere();

private:
    uint32_t x_, y_, z_, w_;
};

inline Xorshift::Xorshift(uint32_t seed) : x_(seed),
                                           y_(362436069),
                                           z_(521288629),
                                           w_(88675123) {}

inline void Xorshift::setSeed(uint32_t seed)
{
    x_ = seed;
}

inline uint32_t Xorshift::next()
{
    uint32_t t = x_ ^ (x_ << 11);
    x_ = y_;
    y_ = z_;
    z_ = w_;
    w_ = (w_ ^ (w_ >> 19)) ^ (t ^ (t >> 8));

    return w_;
}

template <typename T>
inline T Xorshift::next()
{
    return static_cast<T>(next()) / std::numeric_limits<uint32_t>::max();
}

template <typename Vec3T>
inline Vec3T Xorshift::nextOnHemisphere()
{
    using Scalar = typename Vec3T::Scalar;

    Scalar z = next<Scalar>();
    Scalar phi = next<Scalar>() * TAU_F;
    Scalar coeff = std::sqrt(1 - z * z);
    Scalar x = coeff * std::cos(phi);
    Scalar y = coeff * std::sin(phi);

    return Vec3T(x, y, z);
}

inline uint32_t Xorshift::random()
{
    static uint32_t x = 123456789;
    static uint32_t y = 362436069;
    static uint32_t z = 521288629;
    static uint32_t w = 88675123;

    uint32_t t = x ^ (x << 11);
    x = y;
    y = z;
    z = w;
    w = (w ^ (w >> 19)) ^ (t ^ (t >> 8));

    return w;
}

template <typename T>
inline T Xorshift::random()
{
    return static_cast<T>(random()) / std::numeric_limits<uint32_t>::max();
}

template <typename Vec3T>
inline Vec3T Xorshift::randomOnHemisphere()
{
    using Scalar = typename Vec3T::Scalar;

    Scalar z = random<Scalar>();
    Scalar phi = random<Scalar>() * TAU_F;
    Scalar coeff = std::sqrt(1 - z * z);
    Scalar x = coeff * std::cos(phi);
    Scalar y = coeff * std::sin(phi);

    return Vec3T(x, y, z);
}

} // namespace lb

#endif // LIBBSDF_XORSHIFT_H
