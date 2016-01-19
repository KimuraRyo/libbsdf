// =================================================================== //
// Copyright (C) 2014-2015 Kimura Ryo                                  //
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
     * \param equalIntervalAngles If this parameter is true, angles of sample points
     * are equally-spaced intervals.
     */
    CoordinatesBrdf(int         numAngles0,
                    int         numAngles1,
                    int         numAngles2,
                    int         numAngles3,
                    ColorModel  colorModel = RGB_MODEL,
                    int         numWavelengths = 3,
                    bool        equalIntervalAngles = false);

    /*! Constructs a BRDF from lb::Brdf and angle lists. */
    CoordinatesBrdf(const Brdf&     brdf,
                    const Arrayf&   angles0,
                    const Arrayf&   angles1,
                    const Arrayf&   angles2,
                    const Arrayf&   angles3);

    /*! Constructs a BRDF from lb::Brdf and the numbers of angles. Angles are equally-spaced intervals. */
    CoordinatesBrdf(const Brdf& brdf,
                    int         numAngles0,
                    int         numAngles1,
                    int         numAngles2,
                    int         numAngles3);

    /*! Copies and constructs a BRDF. */
    CoordinatesBrdf(const CoordinatesBrdf& brdf);

    virtual ~CoordinatesBrdf();

    /*! Gets the spectrum of the BRDF at incoming and outgoing directions. */
    Spectrum getSpectrum(const Vec3& inDir, const Vec3& outDir) const;

    /*! Gets the value of the BRDF at incoming and outgoing directions and the index of wavelength. */
    float getValue(const Vec3& inDir, const Vec3& outDir, int wavelengthIndex) const;

    /*!
     * Computes incoming and outgoing directions of a Cartesian coordinate system
     * using a set of angle indices.
     */
    void getInOutDirection(int index0, int index1, int index2, int index3,
                           Vec3* inDir, Vec3* outDir) const;

    /*!
     * Converts from four angles to incoming and outgoing directions and
     * assigns them to \a inDir and \a outDir.
     */
    void toXyz(float angle0, float angle1, float angle2, float angle3,
               Vec3* inDir, Vec3* outDir) const;

    /*!
     * Converts from incoming and outgoing directions to four angles and
     * assigns them to \a angle0, \a angle1, \a angle2, and \a angle3.
     */
    void fromXyz(const Vec3& inDir, const Vec3& outDir,
                 float* angle0, float* angle1, float* angle2, float* angle3) const;

    /*!
     * Converts from incoming and outgoing directions to three angles for an isotropic BRDF and
     * assigns them to \a angle0, \a angle2, and \a angle3.
     */
    void fromXyz(const Vec3& inDir, const Vec3& outDir,
                 float* angle0, float* angle2, float* angle3) const;

    std::string getAngle0Name() const; /*!< Gets a name of angle0. */
    std::string getAngle1Name() const; /*!< Gets a name of angle1. */
    std::string getAngle2Name() const; /*!< Gets a name of angle2. */
    std::string getAngle3Name() const; /*!< Gets a name of angle3. */

    /*!
     * Expands minimum angles to 0 and maximum angles to MAX_ANGLE,
     * and constructs the extrapolated sample set.
     */
    bool expandAngles();

    /*! Clamps all angles to minimum and maximum values of each coordinate system. */
    void clampAngles();

protected:
    /*! Sets the angle0 at the index. The array of angles must be sorted in ascending order. */
    void setAngle0(int index, float angle);
    void setAngle1(int index, float angle); /*!< Sets the angle1 at the index. \sa setAngle0() */
    void setAngle2(int index, float angle); /*!< Sets the angle2 at the index. \sa setAngle0() */
    void setAngle3(int index, float angle); /*!< Sets the angle3 at the index. \sa setAngle0() */

private:
    /*! Copy operator is disabled. */
    CoordinatesBrdf& operator=(const CoordinatesBrdf&);

    /*! Initializes angle lists consisting of equal interval angles. */
    void initializeEqualIntervalAngles();

    /*! Initializes spectra using lb::Brdf. */
    void initializeSpectra(const Brdf& brdf);
};

template <typename CoordSysT>
CoordinatesBrdf<CoordSysT>::CoordinatesBrdf(int         numAngles0,
                                            int         numAngles1,
                                            int         numAngles2,
                                            int         numAngles3,
                                            ColorModel  colorModel,
                                            int         numWavelengths,
                                            bool        equalIntervalAngles)
                                            : Brdf(numAngles0,
                                                   numAngles1,
                                                   numAngles2,
                                                   numAngles3,
                                                   colorModel,
                                                   numWavelengths)
{
    if (equalIntervalAngles) {
        initializeEqualIntervalAngles();
    }
}

template <typename CoordSysT>
CoordinatesBrdf<CoordSysT>::CoordinatesBrdf(const Brdf&     brdf,
                                            const Arrayf&   angles0,
                                            const Arrayf&   angles1,
                                            const Arrayf&   angles2,
                                            const Arrayf&   angles3)
                                            : Brdf()
{
    const SampleSet* ss = brdf.getSampleSet();
    samples_ = new SampleSet(angles0.size(), angles1.size(), angles2.size(), angles3.size(),
                             ss->getColorModel(), ss->getNumWavelengths());
    samples_->getAngles0() = angles0;
    samples_->getAngles1() = angles1;
    samples_->getAngles2() = angles2;
    samples_->getAngles3() = angles3;
    samples_->getWavelengths() = ss->getWavelengths();

    samples_->updateAngleAttributes();
    initializeSpectra(brdf);
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
}

template <typename CoordSysT>
CoordinatesBrdf<CoordSysT>::CoordinatesBrdf(const CoordinatesBrdf& brdf) : Brdf(brdf) {}

template <typename CoordSysT>
CoordinatesBrdf<CoordSysT>::~CoordinatesBrdf() {}

template <typename CoordSysT>
Spectrum CoordinatesBrdf<CoordSysT>::getSpectrum(const Vec3& inDir, const Vec3& outDir) const
{
    Spectrum sp;
    Sampler::getSpectrum<CoordSysT, LinearInterpolator>(*samples_, inDir, outDir, &sp);
    return sp;
}

template <typename CoordSysT>
float CoordinatesBrdf<CoordSysT>::getValue(const Vec3& inDir, const Vec3& outDir, int wavelengthIndex) const
{
    return Sampler::getValue<CoordSysT, LinearInterpolator>(*samples_, inDir, outDir, wavelengthIndex);
}

template <typename CoordSysT>
void CoordinatesBrdf<CoordSysT>::getInOutDirection(int index0, int index1, int index2, int index3,
                                                   Vec3* inDir, Vec3* outDir) const
{
    CoordSysT::toXyz(samples_->getAngle0(index0),
                     samples_->getAngle1(index1),
                     samples_->getAngle2(index2),
                     samples_->getAngle3(index3),
                     inDir, outDir);

    inDir->normalize();
    outDir->normalize();
}

template <typename CoordSysT>
void CoordinatesBrdf<CoordSysT>::toXyz(float angle0, float angle1, float angle2, float angle3,
                                       Vec3* inDir, Vec3* outDir) const
{
    CoordSysT::toXyz(angle0, angle1, angle2, angle3, inDir, outDir);
}

template <typename CoordSysT>
void CoordinatesBrdf<CoordSysT>::fromXyz(const Vec3& inDir, const Vec3& outDir,
                                         float* angle0, float* angle1, float* angle2, float* angle3) const
{
    CoordSysT::fromXyz(inDir, outDir, angle0, angle1, angle2, angle3);
}

template <typename CoordSysT>
void CoordinatesBrdf<CoordSysT>::fromXyz(const Vec3& inDir, const Vec3& outDir,
                                         float* angle0, float* angle2, float* angle3) const
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
bool CoordinatesBrdf<CoordSysT>::expandAngles()
{
    bool expanded = false;

    Arrayf angles0 = samples_->getAngles0();
    Arrayf angles1 = samples_->getAngles1();
    Arrayf angles2 = samples_->getAngles2();
    Arrayf angles3 = samples_->getAngles3();

    if (angles0[0] != 0.0f) { appendElement(&angles0, 0.0f); }
    if (angles1[0] != 0.0f) { appendElement(&angles1, 0.0f); }
    if (angles2[0] != 0.0f) { appendElement(&angles2, 0.0f); }
    if (angles3[0] != 0.0f) { appendElement(&angles3, 0.0f); }

    const float maxAngle0 = CoordSysT::MAX_ANGLE0;
    const float maxAngle1 = CoordSysT::MAX_ANGLE1;
    const float maxAngle2 = CoordSysT::MAX_ANGLE2;
    const float maxAngle3 = CoordSysT::MAX_ANGLE3;

    if (angles0.size() > 1 && angles0[angles0.size() - 1] != maxAngle0) { appendElement(&angles0, maxAngle0); }
    if (angles1.size() > 1 && angles1[angles1.size() - 1] != maxAngle1) { appendElement(&angles1, maxAngle1); }
    if (angles2.size() > 1 && angles2[angles2.size() - 1] != maxAngle2) { appendElement(&angles2, maxAngle2); }
    if (angles3.size() > 1 && angles3[angles3.size() - 1] != maxAngle3) { appendElement(&angles3, maxAngle3); }

    bool sizeEqual = (angles0.size() != samples_->getNumAngles0() ||
                      angles1.size() != samples_->getNumAngles1() ||
                      angles2.size() != samples_->getNumAngles2() ||
                      angles3.size() != samples_->getNumAngles3());
    if (sizeEqual) {
        std::sort(angles0.data(), angles0.data() + angles0.size());
        std::sort(angles1.data(), angles1.data() + angles1.size());
        std::sort(angles2.data(), angles2.data() + angles2.size());
        std::sort(angles3.data(), angles3.data() + angles3.size());

        CoordinatesBrdf<CoordSysT> origBrdf(*this);

        samples_->resizeAngles(angles0.size(), angles1.size(), angles2.size(), angles3.size());
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
    Arrayf& angles0 = samples_->getAngles0();
    Arrayf& angles1 = samples_->getAngles1();
    Arrayf& angles2 = samples_->getAngles2();
    Arrayf& angles3 = samples_->getAngles3();

    angles0 = angles0.cwiseMax(0.0);
    angles1 = angles1.cwiseMax(0.0);
    angles2 = angles2.cwiseMax(0.0);
    angles3 = angles3.cwiseMax(0.0);

    angles0 = angles0.cwiseMin(CoordSysT::MAX_ANGLE0);
    angles1 = angles1.cwiseMin(CoordSysT::MAX_ANGLE1);
    angles2 = angles2.cwiseMin(CoordSysT::MAX_ANGLE2);
    angles3 = angles3.cwiseMin(CoordSysT::MAX_ANGLE3);
}

template <typename CoordSysT>
inline void CoordinatesBrdf<CoordSysT>::setAngle0(int index, float angle)
{
    samples_->setAngle0(index, clamp(angle, 0.0f, CoordSysT::MAX_ANGLE0));
}

template <typename CoordSysT>
inline void CoordinatesBrdf<CoordSysT>::setAngle1(int index, float angle)
{
    samples_->setAngle1(index, clamp(angle, 0.0f, CoordSysT::MAX_ANGLE1));
}

template <typename CoordSysT>
inline void CoordinatesBrdf<CoordSysT>::setAngle2(int index, float angle)
{
    samples_->setAngle2(index, clamp(angle, 0.0f, CoordSysT::MAX_ANGLE2));
}

template <typename CoordSysT>
inline void CoordinatesBrdf<CoordSysT>::setAngle3(int index, float angle)
{
    samples_->setAngle3(index, clamp(angle, 0.0f, CoordSysT::MAX_ANGLE3));
}

template <typename CoordSysT>
void CoordinatesBrdf<CoordSysT>::initializeEqualIntervalAngles()
{
    samples_->getAngles0() = Arrayf::LinSpaced(samples_->getNumAngles0(), 0.0, CoordSysT::MAX_ANGLE0);
    samples_->getAngles1() = Arrayf::LinSpaced(samples_->getNumAngles1(), 0.0, CoordSysT::MAX_ANGLE1);
    samples_->getAngles2() = Arrayf::LinSpaced(samples_->getNumAngles2(), 0.0, CoordSysT::MAX_ANGLE2);
    samples_->getAngles3() = Arrayf::LinSpaced(samples_->getNumAngles3(), 0.0, CoordSysT::MAX_ANGLE3);

    if (samples_->getNumAngles0() == 1) { setAngle0(0, 0.0f); }
    if (samples_->getNumAngles1() == 1) { setAngle1(0, 0.0f); }
    if (samples_->getNumAngles2() == 1) { setAngle2(0, 0.0f); }
    if (samples_->getNumAngles3() == 1) { setAngle3(0, 0.0f); }

    samples_->updateAngleAttributes();
}

template <typename CoordSysT>
void CoordinatesBrdf<CoordSysT>::initializeSpectra(const Brdf& brdf)
{
    for (int i0 = 0; i0 < samples_->getNumAngles0(); ++i0) {
    for (int i1 = 0; i1 < samples_->getNumAngles1(); ++i1) {
    for (int i2 = 0; i2 < samples_->getNumAngles2(); ++i2) {
    for (int i3 = 0; i3 < samples_->getNumAngles3(); ++i3) {
        Vec3 inDir, outDir;
        getInOutDirection(i0, i1, i2, i3, &inDir, &outDir);
        fixDownwardDir(&inDir);
        fixDownwardDir(&outDir);

        Spectrum sp = brdf.getSpectrum(inDir, outDir);
        samples_->setSpectrum(i0, i1, i2, i3, sp.cwiseMax(0.0));
    }}}}
}

} // namespace lb

#endif // LIBBSDF_COORDINATES_BRDF_H
