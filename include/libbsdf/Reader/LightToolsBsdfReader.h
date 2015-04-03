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
 * \brief The LightToolsBsdfReader class provides the reader of a LightTools BSDF file.
 *
 * A LightTools BSDF file can contain BRDF and BTDF of front and back.
 */
class LightToolsBsdfReader
{
public:
    /*! Reads a LightTools BSDF file and creates the two-sided material of a spherical coordinate system. */
    static TwoSidedMaterial* read(const std::string& fileName);

private:
    enum SymmetryType {
        ASYMMETRICAL
    };

    enum SideType {
        FRONT,
        BACK
    };

    enum ScatterType {
        BRDF,
        BTDF
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
        ScatterType             scatterType;
        TristimulusValueType    tristimulusValueType;
        
        Arrayf samples;

        static bool cmp(Data* lhs, Data* rhs);
    };

    /*! Skips comment lines. */
    static void ignoreCommentLines(std::ifstream& fin);

    /*! Creates a BRDF from a LightTools BRDF data. */
    static SphericalCoordinatesBrdf* createBrdf(std::vector<Data*>&         brdfData,
                                                const std::vector<float>&   outThetaDegrees,
                                                const std::vector<float>&   outPhiDegrees,
                                                ColorModel::Type            colorModel);

    /*! Rotates a BRDF using an outgoing azimuthal angle. */
    static SphericalCoordinatesBrdf* rotateOutPhi(const SphericalCoordinatesBrdf& brdf, float rotationAngle);
};

inline void LightToolsBsdfReader::ignoreCommentLines(std::ifstream& fin)
{
    reader_utility::ignoreCommentLines(fin, "#");
}

} // namespace lb

#endif // LIBBSDF_LIGHTTOOLS_BSDF_READER_H
