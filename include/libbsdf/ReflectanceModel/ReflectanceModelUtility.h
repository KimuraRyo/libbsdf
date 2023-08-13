// =================================================================== //
// Copyright (C) 2015-2023 Kimura Ryo                                  //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#ifndef LIBBSDF_REFLECTANCE_MODEL_UTILITY_H
#define LIBBSDF_REFLECTANCE_MODEL_UTILITY_H

#include <libbsdf/Brdf/Brdf.h>
#include <libbsdf/ReflectanceModel/ReflectanceModel.h>

namespace lb {

/*!
 * \class   ReflectanceModelUtility
 * \brief   The ReflectanceModelUtility class provides the functions for reflectance models.
 */
class ReflectanceModelUtility
{
public:
    /*! Sets up a lb::Brdf using an analytic reflectance or transmittance model. */
    static bool setupBrdf(const ReflectanceModel& model, Brdf* brdf, DataType dataType = BRDF_DATA);

    /*!
     * Sets up a lb::Brdf using an analytic reflectance or transmittance model.
     * Angles are adaptively inserted up to maximum angles.
     *
     * \param numAngles0    Upper limit of the number of angle0.
     * \param numAngles1    Upper limit of the number of angle1.
     * \param numAngles2    Upper limit of the number of angle2.
     * \param numAngles3    Upper limit of the number of angle3.
     */
    static bool setupBrdf(const ReflectanceModel& model,
                          Brdf*                   brdf,
                          int                     numAngles0,
                          int                     numAngles1,
                          int                     numAngles2,
                          int                     numAngles3,
                          DataType                dataType = BRDF_DATA,
                          double                  ior = 1.0);

    /*! Dumps information about parameters to lbInfo. */
    static void dumpParametersInfo(const ReflectanceModel& model);

private:
    static Spectrum getSpectrum(const ReflectanceModel& model,
                                const Brdf&             brdf,
                                DataType                dataType,
                                double                  angle0,
                                double                  angle1,
                                double                  angle2,
                                double                  angle3);

    static bool insertAngle0(const ReflectanceModel&    model,
                             Brdf*                      brdf,
                             int                        numAngles,
                             DataType                   dataType);

    static bool insertAngle1(const ReflectanceModel&    model,
                             Brdf*                      brdf,
                             int                        numAngles,
                             DataType                   dataType);

    static bool insertAngle2(const ReflectanceModel&    model,
                             Brdf*                      brdf,
                             int                        numAngles,
                             DataType                   dataType);

    static bool insertAngle3(const ReflectanceModel&    model,
                             Brdf*                      brdf,
                             int                        numAngles,
                             DataType                   dataType);
};

} // namespace lb

#endif // LIBBSDF_REFLECTANCE_MODEL_UTILITY_H
