// =================================================================== //
// Copyright (C) 2014-2015 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#ifndef LIBBSDF_ASTM_READER_H
#define LIBBSDF_ASTM_READER_H

#include <map>
#include <string>
#include <vector>

#include <libbsdf/Brdf/SphericalCoordinatesBrdf.h>
#include <libbsdf/Common/Global.h>

namespace lb {

/*!
 * \class AstmReader
 * \brief The AstmReader class provides the reader for an ASTM E1392-96(2002) file.
 *
 * lb::SphericalCoordinatesBrdf is created from loaded data.
 * 
 * File format:
 * http://www.astm.org/Standards/E1392.htm
 */
class AstmReader
{
public:
    /*! Reads an ASTM file and creates the BRDF of a spherical coordinate system. */
    static SphericalCoordinatesBrdf* read(const std::string& fileName);

private:
    typedef std::vector<float> AngleList;
    typedef std::map<AngleList, Spectrum, std::less<AngleList>,
                     Eigen::aligned_allocator<std::pair<AngleList, Spectrum> > > SampleMap;

    /*! Fills omitted data using a plane symmetry. */
    static SphericalCoordinatesBrdf* fillSymmetricBrdf(SphericalCoordinatesBrdf* brdf);

    /*! Finds the nearest sample and returns the spectrum. */
    static const Spectrum& findNearestSample(const SampleMap& samples, const AngleList& sampleAngles);
};

} // namespace lb

#endif // LIBBSDF_ASTM_READER_H
