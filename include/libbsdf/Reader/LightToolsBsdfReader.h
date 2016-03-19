// =================================================================== //
// Copyright (C) 2015-2016 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#ifndef LIBBSDF_LIGHTTOOLS_BSDF_READER_H
#define LIBBSDF_LIGHTTOOLS_BSDF_READER_H

#include <string>

#include <libbsdf/Brdf/SphericalCoordinatesBrdf.h>
#include <libbsdf/Brdf/TwoSidedMaterial.h>
#include <libbsdf/Reader/ReaderUtility.h>

namespace lb {

/*!
 * \class LightToolsBsdfReader
 * \brief The LightToolsBsdfReader class provides the reader for a LightTools BSDF file.
 *
 * A LightTools BSDF file can contain BRDFs and BTDFs of front and back.
 * lb::TwoSidedMaterial is created from loaded data.
 */
class LightToolsBsdfReader
{
public:
    /*! Reads a LightTools BSDF file and creates the two-sided material of a spherical coordinate system. */
    static TwoSidedMaterial* read(const std::string& fileName);

private:
    enum SymmetryType {
        UNKNOWN_SYMMETRY = 0,
        ASYMMETRICAL
    };

    enum SideType {
        FRONT_SIDE,
        BACK_SIDE
    };

    enum TristimulusValueType {
        TRIS_X = 0,
        TRIS_Y = 1,
        TRIS_Z = 2
    };

    /*! The data structure of a LightTools BRDF. */
    struct DataBlock
    {
        DataBlock();

        float aoi;          /*!< Angle Of Incidence. */
        float poi;          /*!< Plane Of Incidence. */
        float wavelength;
        float tis;          /*!< Total Integrated Scatter. */

        SideType                sideType;
        DataType                dataType;
        TristimulusValueType    tristimulusValueType;
        
        Arrayf samples;

        static bool cmp(DataBlock* lhs, DataBlock* rhs);
    };

    /*! Skips comment lines. */
    static void ignoreCommentLines(std::istream& stream);

    /*! Creates a BRDF from LightTools BRDF data. */
    static SphericalCoordinatesBrdf* createBrdf(std::vector<DataBlock*>&    brdfData,
                                                const std::vector<float>&   outThetaDegrees,
                                                const std::vector<float>&   outPhiDegrees,
                                                ColorModel                  colorModel);
};

inline void LightToolsBsdfReader::ignoreCommentLines(std::istream& stream)
{
    reader_utility::ignoreCommentLines(stream, "#");
}

} // namespace lb

#endif // LIBBSDF_LIGHTTOOLS_BSDF_READER_H
