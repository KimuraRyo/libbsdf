// =================================================================== //
// Copyright (C) 2014-2023 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#ifndef LIBBSDF_COORDINATES_BRDF_H
#define LIBBSDF_COORDINATES_BRDF_H

#include <libbsdf/Brdf/Brdf.h>
#include <libbsdf/Brdf/LinearInterpolator.h>
#include <libbsdf/Brdf/Sampler.h>

namespace lb {

/*!
 * \class   CoordinatesBrdf
 * \brief   The CoordinatesBrdf class provides the BRDF using the template of a coordinate system.
 *
 * The properties and data of the BRDF can be get and set with Brdf::samples_.
 */
template <typename CoordSysT>
class CoordinatesBrdf : public Brdf
{
public:
    /*!
     * Constructs a BRDF.
     *
     * \param equalIntervalAngles If this parameter is true, angles of sample points are equally-spaced intervals.
     */
    CoordinatesBrdf(int        numAngles0,
                    int        numAngles1,
                    int        numAngles2,
                    int        numAngles3,
                    ColorModel colorModel = RGB_MODEL,
                    int        numWavelengths = 3,
                    bool       equalIntervalAngles = false);

    /*! Constructs a BRDF from lb::Brdf and angle lists. */
    CoordinatesBrdf(const Brdf&   brdf,
                    const Arrayd& angles0,
                    const Arrayd& angles1,
                    const Arrayd& angles2,
                    const Arrayd& angles3);

    /*! Constructs a BRDF from lb::Brdf and the numbers of angles. Angles are equally-spaced intervals. */
    CoordinatesBrdf(const Brdf& brdf,
                    int         numAngles0,
                    int         numAngles1,
                    int         numAngles2,
                    int         numAngles3);

    /*! Constructs an empty BRDF. Brdf::samples_ must be initialized in a derived class. */
    CoordinatesBrdf();

    /*! Copies and constructs a BRDF. */
    CoordinatesBrdf(const CoordinatesBrdf& brdf);

    virtual ~CoordinatesBrdf();

    /*! Virtual copy constructor. */
    CoordinatesBrdf<CoordSysT>* clone() const override;

    /*! Gets the spectrum of the BRDF at incoming and outgoing directions. */
    Spectrum getSpectrum(const Vec3& inDir, const Vec3& outDir) const override;

    /*! Gets the value of the BRDF at incoming and outgoing directions and the index of wavelength. */
    float getValue(const Vec3& inDir, const Vec3& outDir, int wavelengthIndex) const override;

    /*!
     * Computes incoming and outgoing directions of a Cartesian coordinate system
     * using a set of angle indices.
     */
    void getInOutDirection(int   index0,
                           int   index1,
                           int   index2,
                           int   index3,
                           Vec3* inDir,
                           Vec3* outDir) const override;

    /*!
     * Converts from four angles to incoming and outgoing directions and
     * assigns them to \a inDir and \a outDir.
     */
    void toXyz(double angle0,
               double angle1,
               double angle2,
               double angle3,
               Vec3*  inDir,
               Vec3*  outDir) const override;

    /*!
     * Converts from incoming and outgoing directions to four angles and
     * assigns them to \a angle0, \a angle1, \a angle2, and \a angle3.
     */
    void fromXyz(const Vec3& inDir,
                 const Vec3& outDir,
                 double*     angle0,
                 double*     angle1,
                 double*     angle2,
                 double*     angle3) const override;

    /*!
     * Converts from incoming and outgoing directions to three angles for an isotropic BRDF and
     * assigns them to \a angle0, \a angle2, and \a angle3.
     */
    void fromXyz(const Vec3& inDir,
                 const Vec3& outDir,
                 double*     angle0,
                 double*     angle2,
                 double*     angle3) const override;

    std::string getAngle0Name() const override; /*!< Gets a name of angle0. */
    std::string getAngle1Name() const override; /*!< Gets a name of angle1. */
    std::string getAngle2Name() const override; /*!< Gets a name of angle2. */
    std::string getAngle3Name() const override; /*!< Gets a name of angle3. */

    /*!
     * Validates spectra, angles, wavelengths, and other attributes.
     * False is returned if the data contains one of the following:
     *   - Infinite or NaN spectrum
     *   - Negative spectrum on a visible hemisphere
     *   - Outside, infinite, or NaN angle
     *   - Negative, infinite, or NaN wavelength
     *
     * \param verbose If this parameter is true, all warnings of spectra are output.
     */
    bool validate(bool verbose = false) const override;

    /*!
     * Expands minimum angles to MIN_ANGLE and maximum angles to MAX_ANGLE,
     * and constructs the extrapolated sample set.
     */
    bool expandAngles(bool angle0Expanded = true,
                      bool angle1Expanded = true,
                      bool angle2Expanded = true,
                      bool angle3Expanded = true) override;

    /*! Clamps all angles to minimum and maximum values of each coordinate system. */
    void clampAngles() override;

protected:
    /*! Sets the angle0 at the index. The array of angles must be sorted in ascending order. */
    void setAngle0(int index, double angle);
    void setAngle1(int index, double angle); /*!< Sets the angle1 at the index. \sa setAngle0() */
    void setAngle2(int index, double angle); /*!< Sets the angle2 at the index. \sa setAngle0() */
    void setAngle3(int index, double angle); /*!< Sets the angle3 at the index. \sa setAngle0() */

private:
    /*! Copy operator is disabled. */
    CoordinatesBrdf& operator=(const CoordinatesBrdf&);

    /*! Initializes angle lists consisting of equal interval angles. */
    void initializeEqualIntervalAngles();
};

template <typename CoordSysT>
CoordinatesBrdf<CoordSysT>::CoordinatesBrdf(int        numAngles0,
                                            int        numAngles1,
                                            int        numAngles2,
                                            int        numAngles3,
                                            ColorModel colorModel,
                                            int        numWavelengths,
                                            bool       equalIntervalAngles)
    : Brdf(numAngles0, numAngles1, numAngles2, numAngles3, colorModel, numWavelengths)
{
    if (equalIntervalAngles) {
        initializeEqualIntervalAngles();
    }
}

template <typename CoordSysT>
CoordinatesBrdf<CoordSysT>::CoordinatesBrdf(const Brdf&   brdf,
                                            const Arrayd& angles0,
                                            const Arrayd& angles1,
                                            const Arrayd& angles2,
                                            const Arrayd& angles3)
                                            : Brdf()
{
    const SampleSet* ss = brdf.getSampleSet();
    samples_ = new SampleSet(static_cast<int>(angles0.size()),
                             static_cast<int>(angles1.size()),
                             static_cast<int>(angles2.size()),
                             static_cast<int>(angles3.size()),
                             ss->getColorModel(),
                             ss->getNumWavelengths());
    samples_->getAngles0() = angles0;
    samples_->getAngles1() = angles1;
    samples_->getAngles2() = angles2;
    samples_->getAngles3() = angles3;
    samples_->getWavelengths() = ss->getWavelengths();

    samples_->updateAngleAttributes();
    initializeSpectra(brdf);

    reductionType_ = brdf.getReductionType();
    sourceType_ = brdf.getSourceType();
}

template <typename CoordSysT>
CoordinatesBrdf<CoordSysT>::CoordinatesBrdf(const Brdf& brdf,
                                            int         numAngles0,
                                            int         numAngles1,
                                            int         numAngles2,
                                            int         numAngles3)
    : Brdf()
{
    assert(numAngles0 > 0 && numAngles1 > 0 && numAngles2 > 0 && numAngles3 > 0);

    const SampleSet* ss = brdf.getSampleSet();
    samples_ = new SampleSet(numAngles0, numAngles1, numAngles2, numAngles3,
                             ss->getColorModel(), ss->getNumWavelengths());
    initializeEqualIntervalAngles();
    samples_->getWavelengths() = ss->getWavelengths();

    initializeSpectra(brdf);

    reductionType_ = brdf.getReductionType();
    sourceType_ = brdf.getSourceType();
}

template <typename CoordSysT>
CoordinatesBrdf<CoordSysT>::CoordinatesBrdf() : Brdf() {}

template <typename CoordSysT>
CoordinatesBrdf<CoordSysT>::CoordinatesBrdf(const CoordinatesBrdf& brdf) : Brdf(brdf) {}

template <typename CoordSysT>
CoordinatesBrdf<CoordSysT>::~CoordinatesBrdf() {}

template <typename CoordSysT>
CoordinatesBrdf<CoordSysT>* CoordinatesBrdf<CoordSysT>::clone() const
{
    lbTrace << "[CoordinatesBrdf::clone]";
    return new CoordinatesBrdf<CoordSysT>(*this);
}

template <typename CoordSysT>
Spectrum CoordinatesBrdf<CoordSysT>::getSpectrum(const Vec3& inDir, const Vec3& outDir) const
{
    return Sampler::getSpectrum<CoordSysT, LinearInterpolator>(*samples_, inDir, outDir);
}

template <typename CoordSysT>
float CoordinatesBrdf<CoordSysT>::getValue(const Vec3& inDir, const Vec3& outDir, int wavelengthIndex) const
{
    return Sampler::getValue<CoordSysT, LinearInterpolator>(*samples_, inDir, outDir, wavelengthIndex);
}

template <typename CoordSysT>
void CoordinatesBrdf<CoordSysT>::getInOutDirection(int   index0,
                                                   int   index1,
                                                   int   index2,
                                                   int   index3,
                                                   Vec3* inDir,
                                                   Vec3* outDir) const
{
    toXyz(samples_->getAngle0(index0),
          samples_->getAngle1(index1),
          samples_->getAngle2(index2),
          samples_->getAngle3(index3),
          inDir, outDir);
}

template <typename CoordSysT>
void CoordinatesBrdf<CoordSysT>::toXyz(double angle0,
                                       double angle1,
                                       double angle2,
                                       double angle3,
                                       Vec3*  inDir,
                                       Vec3*  outDir) const
{
    CoordSysT::toXyz(angle0, angle1, angle2, angle3, inDir, outDir);
}

template <typename CoordSysT>
void CoordinatesBrdf<CoordSysT>::fromXyz(const Vec3& inDir,
                                         const Vec3& outDir,
                                         double*     angle0,
                                         double*     angle1,
                                         double*     angle2,
                                         double*     angle3) const
{
    CoordSysT::fromXyz(inDir, outDir, angle0, angle1, angle2, angle3);
}

template <typename CoordSysT>
void CoordinatesBrdf<CoordSysT>::fromXyz(const Vec3& inDir,
                                         const Vec3& outDir,
                                         double*     angle0,
                                         double*     angle2,
                                         double*     angle3) const
{
    CoordSysT::fromXyz(inDir, outDir, angle0, angle2, angle3);
}

template <typename CoordSysT>
std::string CoordinatesBrdf<CoordSysT>::getAngle0Name() const
{
    return CoordSysT::ANGLE0_NAME;
}

template <typename CoordSysT>
std::string CoordinatesBrdf<CoordSysT>::getAngle1Name() const
{
    return CoordSysT::ANGLE1_NAME;
}

template <typename CoordSysT>
std::string CoordinatesBrdf<CoordSysT>::getAngle2Name() const
{
    return CoordSysT::ANGLE2_NAME;
}

template <typename CoordSysT>
std::string CoordinatesBrdf<CoordSysT>::getAngle3Name() const
{
    return CoordSysT::ANGLE3_NAME;
}

template <typename CoordSysT>
bool CoordinatesBrdf<CoordSysT>::validate(bool verbose) const
{
    bool valid = samples_->validate(verbose);
    bool spectraValid = true;

    // Spectra
    for (int i0 = 0; i0 < samples_->getNumAngles0(); ++i0) {
        if (!spectraValid && !verbose) break;
    for (int i1 = 0; i1 < samples_->getNumAngles1(); ++i1) {
        if (!spectraValid && !verbose) break;
    for (int i2 = 0; i2 < samples_->getNumAngles2(); ++i2) {
        if (!spectraValid && !verbose) break;
    for (int i3 = 0; i3 < samples_->getNumAngles3(); ++i3) {
        Vec3 inDir, outDir;
        getInOutDirection(i0, i1, i2, i3, &inDir, &outDir);

        if (outDir.z() < 0) continue;

        const Spectrum& sp = samples_->getSpectrum(i0, i1, i2, i3);

        if (sp.allFinite() && sp.minCoeff() < -EPSILON_F * 10) {
            spectraValid = false;
            lbWarn
                << "[CoordinatesBrdf::validate] The spectrum contains negative value(s) at ("
                << i0 << ", " << i1 << ", " << i2 << ", " << i3 << "):\n\t"
                << sp.format(LB_EIGEN_IO_FMT);

            if (!verbose) break;
        }
    }}}}

    Arrayd angles0 = samples_->getAngles0();
    Arrayd angles1 = samples_->getAngles1();
    Arrayd angles2 = samples_->getAngles2();
    Arrayd angles3 = samples_->getAngles3();

    // Angle arrays
    if (angles0.minCoeff() < CoordSysT::MIN_ANGLE0 ||
        angles0.maxCoeff() > CoordSysT::MAX_ANGLE0 + EPSILON_D * CoordSysT::MAX_ANGLE0) {
        valid = false;
        lbWarn
            << "[CoordinatesBrdf::validate] The angle(s) in angles0 is outside of range:\n\t"
            << angles0.format(LB_EIGEN_IO_FMT);
    }
    else {
        lbInfo << "[CoordinatesBrdf::validate] The array of angle0 is valid.";
    }

    if (angles1.minCoeff() < CoordSysT::MIN_ANGLE1 ||
        angles1.maxCoeff() > CoordSysT::MAX_ANGLE1 + EPSILON_D * CoordSysT::MAX_ANGLE1) {
        valid = false;
        lbWarn
            << "[CoordinatesBrdf::validate] The angle(s) in angles1 is outside of range:\n\t"
            << angles1.format(LB_EIGEN_IO_FMT);
    }
    else {
        lbInfo << "[CoordinatesBrdf::validate] The array of angle1 is valid.";
    }

    if (angles2.minCoeff() < CoordSysT::MIN_ANGLE2 ||
        angles2.maxCoeff() > CoordSysT::MAX_ANGLE2 + EPSILON_D * CoordSysT::MAX_ANGLE2) {
        valid = false;
        lbWarn
            << "[CoordinatesBrdf::validate] The angle(s) in angles2 is outside of range:\n\t"
            << angles2.format(LB_EIGEN_IO_FMT);
    }
    else {
        lbInfo << "[CoordinatesBrdf::validate] The array of angle2 is valid.";
    }

    if (angles3.minCoeff() < CoordSysT::MIN_ANGLE3 ||
        angles3.maxCoeff() > CoordSysT::MAX_ANGLE3 + EPSILON_D * CoordSysT::MAX_ANGLE3) {
        valid = false;
        lbWarn
            << "[CoordinatesBrdf::validate] The angle(s) in angles3 is outside of range:\n\t"
            << angles3.format(LB_EIGEN_IO_FMT);
    }
    else {
        lbInfo << "[CoordinatesBrdf::validate] The array of angle3 is valid.";
    }

    return valid;
}

template <typename CoordSysT>
bool CoordinatesBrdf<CoordSysT>::expandAngles(bool angle0Expanded,
                                              bool angle1Expanded,
                                              bool angle2Expanded,
                                              bool angle3Expanded)
{
    bool expanded = false;

    Arrayd angles0 = samples_->getAngles0();
    Arrayd angles1 = samples_->getAngles1();
    Arrayd angles2 = samples_->getAngles2();
    Arrayd angles3 = samples_->getAngles3();

    using array_util::appendElement;

    if (angle0Expanded && !isEqual(angles0[0], CoordSysT::MIN_ANGLE0)) { appendElement(&angles0, CoordSysT::MIN_ANGLE0); }
    if (angle1Expanded && !isEqual(angles1[0], CoordSysT::MIN_ANGLE1)) { appendElement(&angles1, CoordSysT::MIN_ANGLE1); }
    if (angle2Expanded && !isEqual(angles2[0], CoordSysT::MIN_ANGLE2)) { appendElement(&angles2, CoordSysT::MIN_ANGLE2); }
    if (angle3Expanded && !isEqual(angles3[0], CoordSysT::MIN_ANGLE3)) { appendElement(&angles3, CoordSysT::MIN_ANGLE3); }

    constexpr double maxAngle0 = CoordSysT::MAX_ANGLE0;
    constexpr double maxAngle1 = CoordSysT::MAX_ANGLE1;
    constexpr double maxAngle2 = CoordSysT::MAX_ANGLE2;
    constexpr double maxAngle3 = CoordSysT::MAX_ANGLE3;

    if (angle0Expanded && !isEqual(angles0[angles0.size() - 1], maxAngle0)) { appendElement(&angles0, maxAngle0); }
    if (angle2Expanded && !isEqual(angles2[angles2.size() - 1], maxAngle2)) { appendElement(&angles2, maxAngle2); }
    if (angle3Expanded && !isEqual(angles3[angles3.size() - 1], maxAngle3)) { appendElement(&angles3, maxAngle3); }

    if (angle1Expanded &&
        !samples_->isIsotropic() &&
        !isEqual(angles1[angles1.size() - 1], maxAngle1)) {
        appendElement(&angles1, maxAngle1);
    }

    int numAngles0 = static_cast<int>(angles0.size());
    int numAngles1 = static_cast<int>(angles1.size());
    int numAngles2 = static_cast<int>(angles2.size());
    int numAngles3 = static_cast<int>(angles3.size());

    bool angleAppended = (numAngles0 != samples_->getNumAngles0() ||
                          numAngles1 != samples_->getNumAngles1() ||
                          numAngles2 != samples_->getNumAngles2() ||
                          numAngles3 != samples_->getNumAngles3());
    if (angleAppended) {
        std::sort(angles0.data(), angles0.data() + numAngles0);
        std::sort(angles1.data(), angles1.data() + numAngles1);
        std::sort(angles2.data(), angles2.data() + numAngles2);
        std::sort(angles3.data(), angles3.data() + numAngles3);

        CoordinatesBrdf<CoordSysT> origBrdf(*this);

        samples_->resizeAngles(numAngles0, numAngles1, numAngles2, numAngles3);
        samples_->getAngles0() = angles0;
        samples_->getAngles1() = angles1;
        samples_->getAngles2() = angles2;
        samples_->getAngles3() = angles3;

        samples_->updateAngleAttributes();
        initializeSpectra(origBrdf);

        expanded = true;
    }

    return expanded;
}

template <typename CoordSysT>
void CoordinatesBrdf<CoordSysT>::clampAngles()
{
    Arrayd& angles0 = samples_->getAngles0();
    Arrayd& angles1 = samples_->getAngles1();
    Arrayd& angles2 = samples_->getAngles2();
    Arrayd& angles3 = samples_->getAngles3();

    angles0 = angles0.cwiseMax(CoordSysT::MIN_ANGLE0);
    angles1 = angles1.cwiseMax(CoordSysT::MIN_ANGLE1);
    angles2 = angles2.cwiseMax(CoordSysT::MIN_ANGLE2);
    angles3 = angles3.cwiseMax(CoordSysT::MIN_ANGLE3);

    angles0 = angles0.cwiseMin(CoordSysT::MAX_ANGLE0);
    angles1 = angles1.cwiseMin(CoordSysT::MAX_ANGLE1);
    angles2 = angles2.cwiseMin(CoordSysT::MAX_ANGLE2);
    angles3 = angles3.cwiseMin(CoordSysT::MAX_ANGLE3);
}

template <typename CoordSysT>
void CoordinatesBrdf<CoordSysT>::setAngle0(int index, double angle)
{
    samples_->setAngle0(index, clamp(angle, CoordSysT::MIN_ANGLE0, CoordSysT::MAX_ANGLE0));
}

template <typename CoordSysT>
void CoordinatesBrdf<CoordSysT>::setAngle1(int index, double angle)
{
    samples_->setAngle1(index, clamp(angle, CoordSysT::MIN_ANGLE1, CoordSysT::MAX_ANGLE1));
}

template <typename CoordSysT>
void CoordinatesBrdf<CoordSysT>::setAngle2(int index, double angle)
{
    samples_->setAngle2(index, clamp(angle, CoordSysT::MIN_ANGLE2, CoordSysT::MAX_ANGLE2));
}

template <typename CoordSysT>
void CoordinatesBrdf<CoordSysT>::setAngle3(int index, double angle)
{
    samples_->setAngle3(index, clamp(angle, CoordSysT::MIN_ANGLE3, CoordSysT::MAX_ANGLE3));
}

template <typename CoordSysT>
void CoordinatesBrdf<CoordSysT>::initializeEqualIntervalAngles()
{
    samples_->getAngles0() =
        Arrayd::LinSpaced(samples_->getNumAngles0(), CoordSysT::MIN_ANGLE0, CoordSysT::MAX_ANGLE0);
    samples_->getAngles1() =
        Arrayd::LinSpaced(samples_->getNumAngles1(), CoordSysT::MIN_ANGLE1, CoordSysT::MAX_ANGLE1);
    samples_->getAngles2() =
        Arrayd::LinSpaced(samples_->getNumAngles2(), CoordSysT::MIN_ANGLE2, CoordSysT::MAX_ANGLE2);
    samples_->getAngles3() =
        Arrayd::LinSpaced(samples_->getNumAngles3(), CoordSysT::MIN_ANGLE3, CoordSysT::MAX_ANGLE3);

    if (samples_->getNumAngles0() == 1) { setAngle0(0, 0); }
    if (samples_->getNumAngles1() == 1) { setAngle1(0, 0); }
    if (samples_->getNumAngles2() == 1) { setAngle2(0, 0); }
    if (samples_->getNumAngles3() == 1) { setAngle3(0, 0); }

    samples_->updateAngleAttributes();
}

} // namespace lb

#endif // LIBBSDF_COORDINATES_BRDF_H
