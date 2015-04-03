// =================================================================== //
// Copyright (C) 2014-2015 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#ifndef LIBBSDF_XORSHIFT_H
#define LIBBSDF_XORSHIFT_H

#include <limits>

namespace lb
{
#if defined(__C99__) || (defined(__GNUC__) && __GNUC__ >= 3)
    #include <inttypes.h>
    #include <stdint.h>
#else
    #if defined(__GNUC__)
        typedef short               int16_t;
        typedef unsigned short      uint16_t;
        typedef int                 int32_t;
        typedef unsigned int        uint32_t;
        typedef long long           int64_t;
        typedef unsigned long long  uint64_t;
    #elif defined(_MSC_VER) || defined(__BORLANDC__)
        typedef __int16             int16_t;
        typedef unsigned __int16    uint16_t;
        typedef __int32             int32_t;
        typedef unsigned __int32    uint32_t;
        typedef __int64             int64_t;
        typedef unsigned __int64    uint64_t;
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

    /*! Generates a random integer. Range is [0,std::numeric_limits<uint32_t>::max()]. */
    static uint32_t random();

    /*! Generates a random floating-point number. Range is [0.0,1.0]. */
    template <typename T>
    static T random();

    /*! Generates a random point on the surface of a unit hemisphere. Z-up coordinate system is used. */
    template <typename Vec3T>
    static Vec3T randomOnHemisphere();
    
private:
    uint32_t x_, y_, z_, w_;
};

inline Xorshift::Xorshift(uint32_t seed) : x_(seed),
                                           y_(362436069),
                                           z_(521288629),
                                           w_(88675123) {}

inline void Xorshift::setSeed(uint32_t seed) { x_ = seed; }

inline uint32_t Xorshift::next()
{
    uint32_t t = x_ ^ (x_ << 11);
    x_ = y_; y_ = z_; z_ = w_;
    return w_ = (w_ ^ (w_ >> 19)) ^ (t ^ (t >> 8));
}

inline uint32_t Xorshift::random()
{
    static uint32_t x = 123456789;
    static uint32_t y = 362436069;
    static uint32_t z = 521288629;
    static uint32_t w = 88675123;
    uint32_t t;

    t = x ^ (x << 11);
    x = y; y = z; z = w;
    return w = (w ^ (w >> 19)) ^ (t ^ (t >> 8));
}

template <typename T>
inline T Xorshift::random()
{
    return static_cast<T>(random()) / static_cast<T>(std::numeric_limits<uint32_t>::max());
}

template <typename Vec3T>
inline Vec3T Xorshift::randomOnHemisphere()
{
    float z = Xorshift::random<float>();
    float phi = Xorshift::random<float>() * 2.0f * PI_F;
    float coeff = std::sqrt(1.0f - z * z);

    return Vec3T(coeff * std::cos(phi), coeff * std::sin(phi), z);
}

} // namespace lb

#endif // LIBBSDF_XORSHIFT_H
