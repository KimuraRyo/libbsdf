// =================================================================== //
// Copyright (C) 2014-2015 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

/*!
 * \file    Spectrum.h
 * \brief   Data types for a spectrum.
 */

#ifndef LIBBSDF_SPECTRUM_H
#define LIBBSDF_SPECTRUM_H

#include <vector>

#include <Eigen/Core>

namespace lb {

/*! \brief The data type of a spectrum. */
typedef Eigen::ArrayXf Spectrum;

/*! \brief The data type of spectra. */
typedef std::vector<Spectrum, Eigen::aligned_allocator<Spectrum> > SpectrumList;

/*! \brief The color model of spectra. */
namespace ColorModel {

enum Type {
    MONOCHROME,
    RGB,
    XYZ,
    SPECTRAL
};

} // namespace ColorModel

} // namespace lb

#endif // LIBBSDF_SPECTRUM_H
