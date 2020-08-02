// =================================================================== //
// Copyright (C) 2020 Kimura Ryo                                       //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#ifndef LIBBSDF_SSDD_WRITER_H
#define LIBBSDF_SSDD_WRITER_H

#include <string>

#include <libbsdf/Brdf/Material.h>

namespace lb {

/*!
 * \class SsddWriter
 * \brief The SsddWriter class provides the writer of an SSDD file.
 */
class SsddWriter
{
public:
    /*! Format type of tabular data. */
    enum class DataFormat
    {
        ASCII_DATA,
        BINARY_DATA
    };

    /*! Writes BRDF, BTDF, specular reflectance, and specular transmittance in an SSDD file. */
    static bool write(const std::string&    fileName,
                      const Material&       material,
                      DataFormat            format = DataFormat::ASCII_DATA,
                      const std::string&    comments = "");

    /*! Writes a BRDF in an SSDD file. */
    static bool write(const std::string&    fileName,
                      const Brdf&           brdf,
                      DataFormat            format = DataFormat::ASCII_DATA,
                      const std::string&    comments = "");

    /*! Writes a BTDF in an SSDD file. */
    static bool write(const std::string&    fileName,
                      const Btdf&           btdf,
                      DataFormat            format = DataFormat::ASCII_DATA,
                      const std::string&    comments = "");

    /*! Writes specular reflectance or specular transmittance in an SSDD file. */
    static bool write(const std::string&    fileName,
                      const SampleSet2D&    specularReflectances,
                      DataType              dataType,
                      DataFormat            format = DataFormat::ASCII_DATA,
                      const std::string&    comments = "");

private:
    /*! Outputs the BRDF data block of an SSDD file to a stream. */
    static bool output(const Brdf& brdf, DataFormat format, std::ostream& stream);

    /*! Outputs the specular reflectance data block of an SSDD file to a stream. */
    static bool output(const SampleSet2D& ss2, DataFormat format, std::ostream& stream);

    /*! Outputs the color model of an SSDD file to a stream. */
    static bool output(const ColorModel&    colorModel,
                       const Arrayf&        wavelengths,
                       std::ostream&        stream);

    /*! Outputs BRDF data in ascii form. */
    static void outputAsciiData(const SampleSet&    ss,
                                std::ostream&       stream);

    /*! Outputs BRDF data in binary form. */
    static void outputBinaryData(const SampleSet&   ss,
                                 std::ostream&      stream);

    /*! Outputs specular reflectance data in ascii form. */
    static void outputAsciiData(const SampleSet2D&  ss2,
                                std::ostream&       stream);

    /*! Outputs specular reflectance data in binary form. */
    static void outputBinaryData(const SampleSet2D& ss2,
                                 std::ostream&      stream);
};

} // namespace lb

#endif // LIBBSDF_SSDD_WRITER_H
