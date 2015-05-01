// =================================================================== //
// Copyright (C) 2014-2015 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#include <libbsdf/Brdf/SpecularCoordinatesBrdf.h>

#include <libbsdf/Brdf/Integrator.h>
#include <libbsdf/Brdf/SampleSet2D.h>
#include <libbsdf/Common/PoissonDiskDistributionOnSphere.h>

using namespace lb;

SpecularCoordinatesBrdf::SpecularCoordinatesBrdf(int        numInTheta,
                                                 int        numInPhi,
                                                 int        numSpecTheta,
                                                 int        numSpecPhi,
                                                 ColorModel colorModel,
                                                 int        numWavelengths,
                                                 bool       equalIntervalAngles)
                                                 : BaseBrdf(numInTheta,
                                                            numInPhi,
                                                            numSpecTheta,
                                                            numSpecPhi,
                                                            colorModel,
                                                            numWavelengths,
                                                            equalIntervalAngles) {}

SpecularCoordinatesBrdf::SpecularCoordinatesBrdf(const Brdf&    brdf,
                                                 const Arrayf&  inThetaAngles,
                                                 const Arrayf&  inPhiAngles,
                                                 const Arrayf&  specThetaAngles,
                                                 const Arrayf&  specPhiAngles)
                                                 : BaseBrdf(brdf,
                                                            inThetaAngles,
                                                            inPhiAngles,
                                                            specThetaAngles,
                                                            specPhiAngles) {}

SpecularCoordinatesBrdf::SpecularCoordinatesBrdf(const Brdf&    brdf,
                                                 int            numInTheta,
                                                 int            numInPhi,
                                                 int            numSpecTheta,
                                                 int            numSpecPhi)
                                                 : BaseBrdf(brdf,
                                                            numInTheta,
                                                            numInPhi,
                                                            numSpecTheta,
                                                            numSpecPhi) {}

SpecularCoordinatesBrdf::SpecularCoordinatesBrdf(const SpecularCoordinatesBrdf& brdf)
                                                 : BaseBrdf(brdf) {}

SpecularCoordinatesBrdf::~SpecularCoordinatesBrdf() {}

void SpecularCoordinatesBrdf::fixEnergyConservation()
{
    SampleSet2D reflectances(getNumInTheta(), getNumInPhi(),
                             samples_->getColorModel(), samples_->getNumWavelengths());
    reflectances.getThetaArray()  = samples_->getAngles0();
    reflectances.getPhiArray()    = samples_->getAngles1();
    reflectances.getWavelengths() = samples_->getWavelengths();

    Integrator integrator(PoissonDiskDistributionOnSphere::NUM_SAMPLES_ON_HEMISPHERE, true);

    for (int inThIndex = 0; inThIndex < getNumInTheta(); ++inThIndex) {
    for (int inPhIndex = 0; inPhIndex < getNumInPhi();   ++inPhIndex) {
        Vec3 inDir = SphericalCoordinateSystem::toXyz(getInTheta(inThIndex), getInPhi(inPhIndex));

        Spectrum sp = integrator.computeReflectance(*this, inDir);
        reflectances.setSpectrum(inThIndex, inPhIndex, sp);

        // Fix samples to conserve energy.
        float maxReflectance = sp.maxCoeff();
        if (maxReflectance > 1.0f) {
            for (int i2 = 0; i2 < samples_->getNumAngles2(); ++i2) {
            for (int i3 = 0; i3 < samples_->getNumAngles3(); ++i3) {
                Spectrum& fixedSp = samples_->getSpectrum(inThIndex, inPhIndex, i2, i3);
                const float coeff = 0.999546f; // Reflectance of Lambertian using lb::Integrator.
                fixedSp /= maxReflectance / coeff;
            }}
        }
    }}
}
