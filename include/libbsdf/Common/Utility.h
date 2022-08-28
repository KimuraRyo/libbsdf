// =================================================================== //
// Copyright (C) 2014-2022 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

/*!
 * \file    Utility.h
 * \brief   The Utility.h header file includes the utility functions for libbsdf.
 */

#ifndef LIBBSDF_UTILITY_H
#define LIBBSDF_UTILITY_H

#include <libbsdf/Common/CieData.h>
#include <libbsdf/Common/Global.h>
#include <libbsdf/Common/Log.h>
#include <libbsdf/Common/SphericalCoordinateSystem.h>

namespace lb {

class Brdf;

/*! \brief Clamps a value between a minimum and maximum value. */
template <typename T>
T clamp(T value, T minValue, T maxValue);

/*! \brief Returns true if two values are nearly equal. */
template <typename T>
bool isEqual(T lhs, T rhs);

/*! \brief Returns -1 for a negative number, +1 for a positive number, and 0 for 0. */
template <typename T>
int sign(T val);

/*! \brief Returns a slightly smaller value. */
template<typename T>
constexpr T decrease(T val);

/*! \brief Returns a slightly larger value. */
template<typename T>
constexpr T increase(T val);

/*! \brief Computes linearly-interpolated values. */
template <typename DataT, typename ParameterT>
DataT lerp(const DataT& v0, const DataT& v1, const ParameterT& t);

/*! \brief Interpolates smoothly between two input values with cubic Hermite interpolation. */
template <typename T>
T smoothstep(const T& v0, const T& v1, const T& t);

/*! \brief Interpolates smoothly between two input values with 5th-order Hermite interpolation. */
template <typename T>
T smootherstep(const T& v0, const T& v1, const T& t);

/*! \brief Computes smoothly interpolated values with cubic Hermite interpolation. */
template <typename DataT, typename ParameterT>
DataT hermiteInterpolation3(const DataT& v0, const DataT& v1, const ParameterT& t);

/*! \brief Computes smoothly interpolated values with 5th-order Hermite interpolation. */
template <typename DataT, typename ParameterT>
DataT hermiteInterpolation5(const DataT& v0, const DataT& v1, const ParameterT& t);

/*! \brief Computes a specular direction. */
template <typename Vec3T>
Vec3T reflect(const Vec3T& dir, const Vec3T& normalDir);

/*! \brief Inverts an outgoing direction using a bilateral symmetry. */
template <typename Vec3T>
Vec3T toBilateralSymmetry(const Vec3T& inDir, const Vec3T& outDir);

/*! \brief Converts a value from radians to degrees. */
template <typename T>
T toDegree(const T& radian);

/*! \brief Converts a value from degrees to radians. */
template <typename T>
T toRadian(const T& degree);

/*! \brief Converts an array from radians to degrees. */
template <typename T>
T toDegrees(const T& radians);

/*! \brief Converts an array from degrees to radians. */
template <typename T>
T toRadians(const T& degrees);

/*! \brief Computes a logarithmically scaled value. */
template <typename ValueT, typename BaseT>
ValueT toLogScale(const ValueT& value, const BaseT& base);

/*! \brief Converts a coordinate system. */
template <typename SrcCoordSysT, typename DestCoordSysT>
void convertCoordinateSystem(float  srcAngle0,
                             float  srcAngle1,
                             float  srcAngle2,
                             float  srcAngle3,
                             float* destAngle0,
                             float* destAngle1,
                             float* destAngle2,
                             float* destAngle3);

/*! \brief Returns true if two sample sets have the same color model and wavelengths. */
template <typename T>
bool hasSameColor(const T& ss0, const T& ss1);

/*! \brief Converts tristimulus values using a color space. */
template <typename Vec3T>
Vec3T convertColorSpace(const Vec3T& values, const float* matrix);

/*! \brief Converts from CIE XYZ to sRGB. */
template <typename Vec3T>
Vec3T xyzToSrgb(const Vec3T& xyz);

/*! \brief Converts from sRGB to CIE XYZ. */
template <typename Vec3T>
Vec3T srgbToXyz(const Vec3T& rgb);

/*! \brief Converts from CIE XYZ to Adobe RGB (1998). */
template <typename Vec3T>
Vec3T xyzToAdobeRgb(const Vec3T& xyz);

/*! \brief Converts from Adobe RGB (1998) to CIE XYZ. */
template <typename Vec3T>
Vec3T adobeRgbToXyz(const Vec3T& rgb);

/*! \brief Converts from CIE XYZ to CIE LAB. */
template <typename Vec3T>
Vec3T xyzToLab(const Vec3T& xyz);

/*! \brief Converts from xyY to CIE XYZ. */
template <typename Vec3T>
Vec3T xyyToXyz(const Vec3T& xyy);

/*! \brief Converts from lb::Spectrum to sRGB. */
template <typename Vec3T, typename DataT>
Vec3T toSrgb(const Spectrum& sp, const DataT& ss);

/*! \brief Computes the color difference using CIEDE2000. */
Vec3::Scalar computeCiede2000(const Vec3& lab0, const Vec3& lab1);

/*! \brief Finds nearest Munsell properties. */
Vec3 findMunsellProperties(const Vec3& xyz, std::string* hue, float* value, int* chroma);

/*! \brief Fixes a direction if the Z-component is negative. */
template <typename Vec3T>
void fixDownwardDir(Vec3T* dir);

/*! \brief Returns true if a direction faces the back of a surface. */
bool isDownwardDir(const Vec3& dir);

/*! \brief Returns true if either an incoming or outgoing direction faces the back of a surface. */
bool hasDownwardDir(const Brdf& brdf, int i0, int i1, int i2, int i3);

/*! \brief Returns an enumerator as an integer. */
template <typename EnumT>
typename std::underlying_type<EnumT>::type asInteger(const EnumT& value);

/*! \brief Returns true if both parameters have common enumerator. */
template <typename EnumT>
bool hasSameEnumerator(const EnumT& value0, const EnumT& value1);

/*! \brief Gets the current date in ISO 8601 format (YYYY-MM-DD). */
std::string getDate();

/*
 * Implementation
 */

template <typename T>
T clamp(T value, T minValue, T maxValue)
{
    using std::max;
    using std::min;
    return max(minValue, min(maxValue, value));
}

template <typename T>
int sign(T val)
{
    return (T(0) < val) - (val < T(0));
}

template <typename T>
constexpr T decrease(T val)
{
    return val - std::numeric_limits<T>::epsilon() * val;
}

template <typename T>
constexpr T increase(T val)
{
    return val + std::numeric_limits<T>::epsilon() * val;
}

template <typename T>
bool isEqual(T lhs, T rhs)
{
    using std::abs;
    using std::max;

    T tolerance = std::numeric_limits<T>::epsilon()
                * max(max(abs(lhs), abs(rhs)), T(1))
                * T(2);
    return (abs(lhs - rhs) <= tolerance);
}

template <typename DataT, typename ParameterT>
DataT lerp(const DataT& v0, const DataT& v1, const ParameterT& t)
{
    return v0 + (v1 - v0) * t;
}

template <typename T>
T smoothstep(const T& v0, const T& v1, const T& t)
{
    T coeff = clamp((t - v0) / (v1 - v0), T(0), T(1));
    return coeff * coeff * (T(3) - T(2) * coeff);
}

template <typename T>
T smootherstep(const T& v0, const T& v1, const T& t)
{
    T coeff = clamp((t - v0) / (v1 - v0), T(0), T(1));
    return coeff * coeff * coeff * (coeff * (coeff * T(6) - T(15)) + T(10));
}

template <typename DataT, typename ParameterT>
DataT hermiteInterpolation3(const DataT& v0, const DataT& v1, const ParameterT& t)
{
    ParameterT coeff = smoothstep(ParameterT(0), ParameterT(1), t);
    return lerp(v0, v1, coeff);
}

template <typename DataT, typename ParameterT>
DataT hermiteInterpolation5(const DataT& v0, const DataT& v1, const ParameterT& t)
{
    ParameterT coeff = smootherstep(ParameterT(0), ParameterT(1), t);
    return lerp(v0, v1, coeff);
}

template <typename T>
bool hasSameColor(const T& ss0, const T& ss1)
{
    bool same = true;

    if (ss0.getColorModel() != ss1.getColorModel()) {
        lbInfo
            << "[lb::hasSameColor] Color models do not match: "
            << ss0.getColorModel() << ", " << ss1.getColorModel();
        same = false;
    }

    if (ss0.getNumWavelengths() != ss1.getNumWavelengths() ||
        !ss0.getWavelengths().isApprox(ss1.getWavelengths())) {
        lbInfo
            << "[lb::hasSameColor] Wavelengths do not match: "
            << ss0.getWavelengths() << ", " << ss1.getWavelengths();
        same = false;
    }

    return same;
}

template <typename Vec3T>
Vec3T reflect(const Vec3T& dir, const Vec3T& normalDir)
{
    return 2 * normalDir.dot(dir) * normalDir - dir;
}

template <typename Vec3T>
Vec3T toBilateralSymmetry(const Vec3T& inDir, const Vec3T& outDir)
{
    using ScalarType = typename Vec3T::Scalar;

    ScalarType outTheta, outPhi;
    SphericalCoordinateSystem::fromXyz(outDir, &outTheta, &outPhi);

    // Compute the outgoing azimuthal angle inverted along the incident plane.
    ScalarType inPhi = SphericalCoordinateSystem::toPhi(inDir);
    ScalarType invertedOutPhi = -outPhi + 2 * inPhi;

    if (invertedOutPhi < ScalarType(0)) {
        invertedOutPhi += ScalarType(TAU_D);
    }
    else if (invertedOutPhi < TAU_D) {
        invertedOutPhi -= ScalarType(TAU_D);
    }

    return SphericalCoordinateSystem::toXyz(outTheta, invertedOutPhi);
}

template <typename T>
T toDegree(const T& radian)
{
    return static_cast<T>(radian * 180.0 / PI_D);
}

template <typename T>
T toRadian(const T& degree)
{
    return static_cast<T>(degree * PI_D / 180.0);
}

template <typename T>
T toDegrees(const T& radians)
{
    using ScalarType = typename T::Scalar;
    return radians * ScalarType(180 / PI_D);
}

template <typename T>
T toRadians(const T& degrees)
{
    using ScalarType = typename T::Scalar;
    return degrees * ScalarType(PI_D / 180);
}

template <typename ValueT, typename BaseT>
ValueT toLogScale(const ValueT& value, const BaseT& base)
{
    using std::log;
    return static_cast<ValueT>(log(value + ValueT(1)) / log(base));
}

template <typename SrcCoordSysT, typename DestCoordSysT>
void convertCoordinateSystem(float  srcAngle0,
                             float  srcAngle1,
                             float  srcAngle2,
                             float  srcAngle3,
                             float* destAngle0,
                             float* destAngle1,
                             float* destAngle2,
                             float* destAngle3)
{
    Vec3 inDir, outDir;
    SrcCoordSysT::toXyz(srcAngle0, srcAngle1, srcAngle2, srcAngle3,
                        &inDir, &outDir);

    DestCoordSysT::fromXyz(inDir, outDir,
                           destAngle0, destAngle1, destAngle2, destAngle3);
}

template <typename Vec3T>
Vec3T convertColorSpace(const Vec3T& values, const float* matrix)
{
    using ScalarType = typename Vec3T::Scalar;

    Eigen::Matrix<ScalarType, 3, 3> mat;
    mat << matrix[0], matrix[1], matrix[2],
           matrix[3], matrix[4], matrix[5],
           matrix[6], matrix[7], matrix[8];

    return mat * values;
}

template <typename Vec3T>
Vec3T xyzToSrgb(const Vec3T& xyz)
{
    return convertColorSpace(xyz, CieData::XYZ_sRGB);
}

template <typename Vec3T>
Vec3T srgbToXyz(const Vec3T& rgb)
{
    return convertColorSpace(rgb, CieData::sRGB_XYZ);
}

template <typename Vec3T>
Vec3T xyzToAdobeRgb(const Vec3T& xyz)
{
    return convertColorSpace(xyz, CieData::XYZ_AdobeRGB);
}

template <typename Vec3T>
Vec3T adobeRgbToXyz(const Vec3T& rgb)
{
    return convertColorSpace(rgb, CieData::AdobeRGB_XYZ);
}

template <typename Vec3T>
Vec3T xyzToLab(const Vec3T& xyz)
{
    using std::pow;

    auto f = [](double t) {
        constexpr double delta = 6.0 / 29.0;
        if (t > pow(delta, 3.0)) {
            return pow(t, 1.0 / 3.0);
        }
        else {
            return t / (3.0 * delta * delta) + (4.0 / 29.0);
        }
    };

    // White point under Illuminant D65
    // https://en.wikipedia.org/wiki/CIELAB_color_space
    constexpr double Xn = 0.950489;
    constexpr double Yn = 1.0;
    constexpr double Zn = 1.088840;

    double fx = f(xyz.x() / Xn);
    double fy = f(xyz.y() / Yn);
    double fz = f(xyz.z() / Zn);

    using ScalarType = typename Vec3T::Scalar;

    ScalarType L = static_cast<ScalarType>(116.0 * fy - 16.0);
    ScalarType a = static_cast<ScalarType>(500.0 * (fx - fy));
    ScalarType b = static_cast<ScalarType>(200.0 * (fy - fz));

    return Vec3T(L, a, b);
}

template <typename Vec3T>
Vec3T xyyToXyz(const Vec3T& xyy)
{
    using ScalarType = typename Vec3T::Scalar;

    ScalarType x = xyy[0];
    ScalarType y = xyy[1];
    ScalarType Y = xyy[2];

    if (y == ScalarType(0)) {
        Vec3T::Zero();
    }

    ScalarType X = Y / y * x;
    ScalarType Z = Y / y * (ScalarType(1) - x - y);

    return Vec3T(X, Y, Z);
}

template <typename Vec3T, typename DataT>
Vec3T toSrgb(const Spectrum& sp, const DataT& ss)
{
    using ScalarType = typename Vec3T::Scalar;

    if (ss.getColorModel() == RGB_MODEL) {
        return sp.cast<ScalarType>();
    }
    else if (ss.getColorModel() == XYZ_MODEL) {
        return xyzToSrgb<Vec3T>(sp.cast<ScalarType>());
    }
    else if (ss.getNumWavelengths() == 1) {
        return Vec3T(sp[0], sp[0], sp[0]);
    }
    else {
        Vec3 rgb = SpectrumUtility::spectrumToSrgb(sp, ss.getWavelengths());
        return rgb.cast<ScalarType>();
    }
}

template <typename Vec3T>
void fixDownwardDir(Vec3T* dir)
{
    Vec3T& d = *dir;
    if (d[2] < 0) {
        d[2] = 0;
        if (d[0] == 0 && d[1] == 0) {
            d[0] = 1;
        }
        else {
            d.normalize();
        }
    }
}

inline bool isDownwardDir(const Vec3& dir)
{
    return (dir.z() < -0.00001);
}

template <typename EnumT>
typename std::underlying_type<EnumT>::type asInteger(const EnumT& value)
{
    return static_cast<typename std::underlying_type<EnumT>::type>(value);
}

template <typename EnumT>
bool hasSameEnumerator(const EnumT& value0, const EnumT& value1)
{
    return static_cast<bool>(asInteger(value0) & asInteger(value1));
}

} // namespace lb

#endif // LIBBSDF_UTILITY_H
