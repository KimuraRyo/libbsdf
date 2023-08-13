// =================================================================== //
// Copyright (C) 2020-2023 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#include <libbsdf/Writer/SsddWriter.h>

#include <fstream>
#include <sstream>

#include <libbsdf/Brdf/HalfDifferenceCoordinatesBrdf.h>
#include <libbsdf/Brdf/SpecularCoordinatesBrdf.h>
#include <libbsdf/Brdf/SphericalCoordinatesBrdf.h>

#include <libbsdf/Common/Utility.h>
#include <libbsdf/Common/Version.h>

#include <libbsdf/Reader/ReaderUtility.h>
#include <libbsdf/Reader/SsddUtility.h>

using namespace lb;

using std::endl;

bool SsddWriter::write(const std::string& fileName,
                       const Material&    material,
                       DataFormat         format,
                       const std::string& comments)
{
    if (!material.getBsdf() &&
        !material.getSpecularReflectances() &&
        !material.getSpecularTransmittances()) {
        lbError << "[SsddWriter::write] Data is empty.";
        return false;
    }

    std::ofstream ofs(fileName.c_str(), std::ios_base::binary);
    if (ofs.fail()) {
        lbError << "[SsddWriter::write] Could not open: " << fileName;
        return false;
    }

    std::ios_base::sync_with_stdio(false);

    if (!comments.empty()) {
        std::stringstream ss(comments);
        std::string lineStr;
        while (!reader_utility::safeGetline(ss, lineStr).eof()) {
            ofs << "# " << lineStr << endl;
        }
    }

    // Output the file header.
    ofs << ssdd::VERSION << " " << ssdd::FILE_VERSION << endl;
    ofs << endl;
    //ofs << ssdd::SOFTWARE << "" << endl;
    ofs << ssdd::API << " libbsdf-" << getVersion() << endl;
    ofs << ssdd::DATE << " " << getDate() << endl;

    std::shared_ptr<const Bsdf> bsdf = material.getBsdf();
    std::shared_ptr<const Brdf> brdf = bsdf->getBrdf();
    std::shared_ptr<const Btdf> btdf = bsdf->getBtdf();
    std::shared_ptr<const SampleSet2D> specR = material.getSpecularReflectances();
    std::shared_ptr<const SampleSet2D> specT = material.getSpecularTransmittances();

    // Output BRDF data.
    if (brdf) {
        ofs << endl << ssdd::DATA_TYPE << " " << ssdd::DATA_TYPE_BRDF << endl;

        if (!output(*brdf, format, ofs)) {
            lbError << "[SsddWriter::write] Failed to save: " << fileName;
            return false;
        }
    }

    // Output BTDF data.
    if (btdf) {
        ofs << endl << ssdd::DATA_TYPE << " " << ssdd::DATA_TYPE_BTDF << endl;

        if (!output(*btdf->getBrdf(), format, ofs)) {
            lbError << "[SsddWriter::write] Failed to save: " << fileName;
            return false;
        }
    }

    // Output specular reflectance data.
    if (specR) {
        ofs << endl << ssdd::DATA_TYPE << " " << ssdd::DATA_TYPE_SPECULAR_REFLECTANCE << endl;

        if (!output(*specR, format, ofs)) {
            lbError << "[SsddWriter::write] Failed to save: " << fileName;
            return false;
        }
    }

    // Output specular transmittance data.
    if (specT) {
        ofs << endl << ssdd::DATA_TYPE << " " << ssdd::DATA_TYPE_SPECULAR_TRANSMITTANCE << endl;

        if (!output(*specT, format, ofs)) {
            lbError << "[SsddWriter::write] Failed to save: " << fileName;
            return false;
        }
    }

    return true;
}

bool SsddWriter::write(const std::string& fileName,
                       const Brdf&        brdf,
                       DataFormat         format,
                       const std::string& comments)
{
    std::shared_ptr<Bsdf> bsdf = std::make_shared<Bsdf>(std::shared_ptr<Brdf>(brdf.clone()), nullptr);
    std::unique_ptr<Material> material(new Material(bsdf));
    return SsddWriter::write(fileName, *material, format, comments);
}

bool SsddWriter::write(const std::string& fileName,
                       const Btdf&        btdf,
                       DataFormat         format,
                       const std::string& comments)
{
    std::shared_ptr<Bsdf> bsdf = std::make_shared<Bsdf>(nullptr, std::make_shared<Btdf>(btdf));
    std::unique_ptr<Material> material(new Material(bsdf));
    return SsddWriter::write(fileName, *material, format, comments);
}

bool SsddWriter::write(const std::string& fileName,
                       const SampleSet2D& specularReflectances,
                       DataType           dataType,
                       DataFormat         format,
                       const std::string& comments)
{
    std::unique_ptr<Material> material;

    switch (dataType) {
        case SPECULAR_REFLECTANCE_DATA:
            material.reset(new Material(nullptr,
                                        std::make_shared<SampleSet2D>(specularReflectances),
                                        nullptr));
            break;
        case SPECULAR_TRANSMITTANCE_DATA:
            material.reset(new Material(nullptr,
                                        nullptr,
                                        std::make_shared<SampleSet2D>(specularReflectances)));
            break;
        default:
            lbError << "[SsddWriter::write] Invalid data type: " << dataType;
            return false;
    }

    return SsddWriter::write(fileName, *material, format, comments);
}

bool SsddWriter::output(const Brdf& brdf, DataFormat format, std::ostream& stream)
{
    const SampleSet* ss = brdf.getSampleSet();

    output(ss->getColorModel(), ss->getWavelengths(), stream);

    stream << ssdd::PARAM_TYPE << " ";
    auto halfDiffBrdf   = dynamic_cast<const HalfDifferenceCoordinatesBrdf*>(&brdf);
    auto specBrdf       = dynamic_cast<const SpecularCoordinatesBrdf*>(&brdf);
    auto spheBrdf       = dynamic_cast<const SphericalCoordinatesBrdf*>(&brdf);
    if (halfDiffBrdf) {
        stream << ssdd::PARAM_TYPE_HALF_DIFF << endl;
    }
    else if (specBrdf) {
        stream << ssdd::PARAM_TYPE_SPECULAR << endl;
    }
    else if (spheBrdf) {
        stream << ssdd::PARAM_TYPE_SPHERICAL << endl;
    }
    else {
        lbError << "[SsddWriter::write] Unknown parameterization type.";
        return false;
    }

    ReductionType reductionType = brdf.getReductionType();
    if (reductionType != ReductionType::NONE) {
        stream << ssdd::REDUCTION_TYPE;

        if (hasSameEnumerator(reductionType, ReductionType::BILATERAL_SYMMETRY)) {
            stream << " " << ssdd::REDUCTION_TYPE_BILATERAL_SYMMETRY;
        }

        if (hasSameEnumerator(reductionType, ReductionType::RECIPROCITY)) {
            stream << " " << ssdd::REDUCTION_TYPE_RECIPROCITY;
        }

        stream << endl;
    }

    Arrayd degrees0 = toDegrees(ss->getAngles0());
    Arrayd degrees1 = toDegrees(ss->getAngles1());
    Arrayd degrees2 = toDegrees(ss->getAngles2());
    Arrayd degrees3 = toDegrees(ss->getAngles3());

    stream << ssdd::PARAM0_LIST << " " << degrees0.format(ssdd::LIST_FORMAT) << endl;
    if (degrees1.size() > 1 || degrees1[0] != 0) {
        stream << ssdd::PARAM1_LIST << " " << degrees1.format(ssdd::LIST_FORMAT) << endl;
    }
    stream << ssdd::PARAM2_LIST << " " << degrees2.format(ssdd::LIST_FORMAT) << endl;
    stream << ssdd::PARAM3_LIST << " " << degrees3.format(ssdd::LIST_FORMAT) << endl;

    if (specBrdf &&
        specBrdf->getNumSpecularOffsets() == specBrdf->getNumInTheta()) {
        Arrayd offsets = toDegrees(specBrdf->getSpecularOffsets());
        stream << ssdd::PARAM4_LIST << " " << offsets.format(ssdd::LIST_FORMAT) << endl;
    }

    switch (format) {
        case DataFormat::ASCII_DATA:
            outputAsciiData(*ss, stream);
            break;
        case DataFormat::BINARY_DATA:
            outputBinaryData(*ss, stream);
            break;
        default:
            lbError << "[SsddWriter::output] Unknown data format: " << asInteger(format);
            return false;
    }

    return true;
}

bool SsddWriter::output(const SampleSet2D& ss2, DataFormat format, std::ostream& stream)
{
    output(ss2.getColorModel(), ss2.getWavelengths(), stream);

    Arrayd degrees0 = toDegrees(ss2.getThetaArray());
    Arrayd degrees1 = toDegrees(ss2.getPhiArray());

    stream << ssdd::PARAM0_LIST << " " << degrees0.format(ssdd::LIST_FORMAT) << endl;
    if (degrees1.size() > 1 || degrees1[0] != 0) {
        stream << ssdd::PARAM1_LIST << " " << degrees1.format(ssdd::LIST_FORMAT) << endl;
    }

    switch (format) {
        case DataFormat::ASCII_DATA:
            outputAsciiData(ss2, stream);
            break;
        case DataFormat::BINARY_DATA:
            outputBinaryData(ss2, stream);
            break;
        default:
            lbError << "[SsddWriter::output] Unknown data format: " << asInteger(format);
            return false;
    }

    return true;
}

bool SsddWriter::output(const ColorModel& colorModel,
                        const Arrayf&     wavelengths,
                        std::ostream&     stream)
{
    stream << ssdd::COLOR_MODEL << " ";

    switch (colorModel) {
        case MONOCHROMATIC_MODEL:
            stream << ssdd::COLOR_MODEL_MONOCHROME << endl;
            break;
        case RGB_MODEL:
            stream << ssdd::COLOR_MODEL_RGB << endl;
            break;
        case XYZ_MODEL:
            stream << ssdd::COLOR_MODEL_XYZ << endl;
            break;
        case SPECTRAL_MODEL:
            stream << ssdd::COLOR_MODEL_SPECTRUM << endl;
            stream << ssdd::WAVELENGTH_LIST << " " << wavelengths.format(ssdd::LIST_FORMAT) << endl;
            break;
        default:
            lbError << "Unknown color model: " << colorModel;
            return false;
    }

    return true;
}

void SsddWriter::outputAsciiData(const SampleSet& ss, std::ostream& stream)
{
    stream << ssdd::DATA << " " << ssdd::DATA_ASCII << endl;

    // Output tabular data.
    for (int i3 = 0; i3 < ss.getNumAngles3(); ++i3) {
        if (ss.getNumAngles3() > 1) {
            stream << "# PARAM3: " << toDegree(ss.getAngle3(i3)) << endl;
        }

        for (int i2 = 0; i2 < ss.getNumAngles2(); ++i2) {
            if (ss.getNumAngles2() > 1) {
                stream << "# PARAM2: " << toDegree(ss.getAngle2(i2)) << endl;
            }

            for (int i1 = 0; i1 < ss.getNumAngles1(); ++i1) {
                if (ss.getNumAngles1() > 1) {
                    stream << "# PARAM1: " << toDegree(ss.getAngle1(i1)) << endl;
                }

                for (int i0 = 0; i0 < ss.getNumAngles0(); ++i0) {
                    stream << ss.getSpectrum(i0, i1, i2, i3).format(ssdd::LIST_FORMAT) << endl;
                }
            }
        }
    }
}

void SsddWriter::outputBinaryData(const SampleSet& ss, std::ostream& stream)
{
    stream << ssdd::DATA << " " << ssdd::DATA_BINARY << endl;

    for (int i3 = 0; i3 < ss.getNumAngles3(); ++i3) {
    for (int i2 = 0; i2 < ss.getNumAngles2(); ++i2) {
    for (int i1 = 0; i1 < ss.getNumAngles1(); ++i1) {
    for (int i0 = 0; i0 < ss.getNumAngles0(); ++i0) {
        // to do: Perform a byte swap in the big-endian system.
        const Spectrum& sp = ss.getSpectrum(i0, i1, i2, i3);
        stream.write(reinterpret_cast<const char*>(sp.data()), sp.size() * sizeof(typename Spectrum::Scalar));
    }}}}
}

void SsddWriter::outputAsciiData(const SampleSet2D& ss2, std::ostream& stream)
{
    stream << ssdd::DATA << " " << ssdd::DATA_ASCII << endl;

    // Output tabular data.
    for (int i1 = 0; i1 < ss2.getNumPhi(); ++i1) {
        if (ss2.getNumPhi() > 1) {
            stream << "# PARAM1: " << toDegree(ss2.getPhi(i1)) << endl;
        }

        for (int i0 = 0; i0 < ss2.getNumTheta(); ++i0) {
            stream << ss2.getSpectrum(i0, i1).format(ssdd::LIST_FORMAT) << endl;
        }
    }
}

void SsddWriter::outputBinaryData(const SampleSet2D& ss2, std::ostream& stream)
{
    stream << ssdd::DATA << " " << ssdd::DATA_BINARY << endl;

    for (int i1 = 0; i1 < ss2.getNumPhi();   ++i1) {
    for (int i0 = 0; i0 < ss2.getNumTheta(); ++i0) {
        // to do: Perform a byte swap in the big-endian system.
        const Spectrum& sp = ss2.getSpectrum(i0, i1);
        stream.write(reinterpret_cast<const char*>(sp.data()), sp.size() * sizeof(typename Spectrum::Scalar));
    }}
}
