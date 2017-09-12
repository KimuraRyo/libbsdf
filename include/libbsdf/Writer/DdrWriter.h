// =================================================================== //
// Copyright (C) 2014-2017 Kimura Ryo                                  //
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
    static bool write(const std::string& fileName, const SpecularCoordinatesBrdf& brdf);

    /*! Converts and exports a DDR or DDT file. */
    static void write(const std::string&    fileName,
                      const Brdf&           brdf,
                      bool                  inDirDependentCoordSysUsed,
                      DataType              dataType);

    /*! Outputs character data of a DDR or DDT file to a stream. */
    static bool output(const SpecularCoordinatesBrdf& brdf, std::ostream& stream);
};

} // namespace lb

#endif // LIBBSDF_DDR_WRITER_H
