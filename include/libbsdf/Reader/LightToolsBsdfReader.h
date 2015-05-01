// =================================================================== //
// Copyright (C) 2015 Kimura Ryo                                       //
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
    struct Data
    {
        Data();

        float aoi;          /*!< Angle Of Incidence. */
        float poi;          /*!< Plane Of Incidence. */
        float wavelength;
        float tis;          /*!< Total Integrated Scatter. */

        SideType                sideType;
        DataType                dataType;
        TristimulusValueType    tristimulusValueType;
        
        Arrayf samples;

        static bool cmp(Data* lhs, Data* rhs);
    };

    /*! Skips comment lines. */
    static void ignoreCommentLines(std::ifstream& fin);

    /*! Creates a BRDF from LightTools BRDF data. */
    static SphericalCoordinatesBrdf* createBrdf(std::vector<Data*>&         brdfData,
                                                const std::vector<float>&   outThetaDegrees,
                                                const std::vector<float>&   outPhiDegrees,
                                                ColorModel                  colorModel);

    /*! Rotates a BRDF using an outgoing azimuthal angle. */
    static SphericalCoordinatesBrdf* rotateOutPhi(const SphericalCoordinatesBrdf& brdf, float rotationAngle);
};

inline void LightToolsBsdfReader::ignoreCommentLines(std::ifstream& fin)
{
    reader_utility::ignoreCommentLines(fin, "#");
}

} // namespace lb

#endif // LIBBSDF_LIGHTTOOLS_BSDF_READER_H
