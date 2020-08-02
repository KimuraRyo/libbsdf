// =================================================================== //
// Copyright (C) 2020 Kimura Ryo                                       //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#ifndef LIBBSDF_SSDD_READER_H
#define LIBBSDF_SSDD_READER_H

#include <string>

#include <libbsdf/Brdf/Material.h>

namespace lb {

/*!
 * \class SsddReader
 * \brief The SsddReader class provides the reader for an SSDD file.
 */
class SsddReader
{
public:
    /*! Reads an SSDD file and creates the material including a BRDF, BTDF, specular reflectance, and specular transmittance. */
    static Material* read(const std::string& fileName);

private:
    struct FileInfo
    {
        std::string version;
        std::string software;
        std::string api;
        std::string date;
    };

    struct DataInfo
    {
        DataInfo() : reductionType(ReductionType::NONE),
                     sourceType(SourceType::UNKNOWN_SOURCE) {};

        DataType            dataType;
        ColorModel          colorModel;
        std::vector<float>  wavelengths;
        std::string         paramType;
        ReductionType       reductionType;
        std::vector<float>  params0, params1, params2, params3, params4;
        std::string         name;
        SourceType          sourceType;
        std::string         device;
        std::string         creationDate;
        std::string         measurementDate;
        std::string         dataMode;
    };

    /*! Reads and creates a BRDF using information. */
    static std::shared_ptr<Brdf> readBrdf(std::ifstream& ifs, const DataInfo& dataInfo);

    /*! Reads and creates specular reflectances using information. */
    static std::shared_ptr<SampleSet2D> readSpecularReflectances(std::ifstream& ifs, const DataInfo& dataInfo);

    /*! Reads BRDF data in ascii form. */
    static bool readAsciiData(std::ifstream& ifs, SampleSet* ss);

    /*! Reads BRDF data in binary form. */
    static bool readBinaryData(std::ifstream& ifs, SampleSet* ss);

    /*! Reads specular reflectance data in ascii form. */
    static bool readAsciiData(std::ifstream& ifs, SampleSet2D* ss2);

    /*! Reads specular reflectance data in binary form. */
    static bool readBinaryData(std::ifstream& ifs, SampleSet2D* ss2);

    /*! Gets float values in a string stream. */
    static std::vector<float> getList(std::stringstream& stream);
};

} // namespace lb

#endif // LIBBSDF_SSDD_READER_H
