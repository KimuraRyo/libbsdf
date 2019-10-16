// =================================================================== //
// Copyright (C) 2014-2019 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#ifndef LIBBSDF_DDR_WRITER_H
#define LIBBSDF_DDR_WRITER_H

#include <string>

#include <libbsdf/Brdf/SpecularCoordinatesBrdf.h>

namespace lb {

/*!
 * \class DdrWriter
 * \brief The DdrWriter class provides the writer of a DDR and DDT file.
 */
class DdrWriter
{
public:
    /*! Writes the BRDF of a specular coordinate system in a DDR or DDT file. */
    static bool write(const std::string&                fileName,
                      const SpecularCoordinatesBrdf&    brdf,
                      const std::string&                comments = "");

    /*!
     * Converts and writes a DDR or DDT file.
     * Samples are extrapolated using \a dataType if it is possible.
     */
    static bool write(const std::string&    fileName,
                      const Brdf&           brdf,
                      DataType              dataType,
                      const std::string&    comments = "");

    /*! Outputs character data of a DDR or DDT file to a stream. */
    static bool output(const SpecularCoordinatesBrdf&   brdf,
                       std::ostream&                    stream,
                       const std::string&               comments = "");

    /*! Converts lb::Brdf to lb::SpecularCoordinatesBrdf. */
    static SpecularCoordinatesBrdf* convert(const Brdf& brdf);

    /*!
     * Arranges a BRDF.
     * Minimum and maximum angles are added, and spectra are extrapolated.
     * BRDF is fixed if the sum of reflectances and transmittances exceed one.
     * Spectra of samples are filled if incoming polar angle is 90 degrees
     */
    static SpecularCoordinatesBrdf* arrange(const SpecularCoordinatesBrdf&  brdf,
                                            DataType                        dataType);
};

} // namespace lb

#endif // LIBBSDF_DDR_WRITER_H
