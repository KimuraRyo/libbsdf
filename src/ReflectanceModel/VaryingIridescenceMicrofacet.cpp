/*
 * The implementation of "A Practical Extension to Microfacet Theory for the
 * Modeling of Varying Iridescence" by Laurent Belcour and Pascal Barla.
 * It is ported from the reference GLSL shader of the supplemental material,
 * which is available at:
 * https://belcour.github.io/blog/research/publication/2017/05/01/brdf-thin-film.html
 */

#include <libbsdf/ReflectanceModel/VaryingIridescenceMicrofacet.h>

#include <cmath>

#include <libbsdf/Common/Utility.h>
#include <libbsdf/ReflectanceModel/Ggx.h>

using namespace lb;

namespace {

double sqr(double x) { return x * x; }

/* Fresnel reflectance and phase shift for the s- and p-polarization of an interface. */
struct FresnelResult
{
    double Rs, Rp;
    double phiS, phiP;
};

/* Fresnel equations for a dielectric/dielectric interface. */
FresnelResult fresnelDielectric(double ct1, double n1, double n2)
{
    double st1 = 1.0 - ct1 * ct1;
    double sqNr = sqr(n1 / n2);

    FresnelResult result;
    if (sqNr * st1 > 1.0) {
        // Total internal reflection.
        double term = std::sqrt(st1 - 1.0 / sqNr) / ct1;
        result.Rs = 1.0;
        result.Rp = 1.0;
        result.phiS = 2.0 * std::atan(-term);
        result.phiP = 2.0 * std::atan(-sqNr * term);
    }
    else {
        double ct2 = std::sqrt(1.0 - sqNr * st1);
        double rs = (n1 * ct1 - n2 * ct2) / (n1 * ct1 + n2 * ct2);
        double rp = (n2 * ct1 - n1 * ct2) / (n2 * ct1 + n1 * ct2);
        result.Rs = rs * rs;
        result.Rp = rp * rp;
        result.phiS = (rs < 0.0) ? PI_D : 0.0;
        result.phiP = (rp < 0.0) ? PI_D : 0.0;
    }

    return result;
}

/* Fresnel equations for a dielectric/conductor interface. */
FresnelResult fresnelConductor(double ct1, double n1, double n2, double k)
{
    if (k == 0.0) {
        return fresnelDielectric(ct1, n1, n2);
    }

    double A = sqr(n2) * (1.0 - sqr(k)) - sqr(n1) * (1.0 - sqr(ct1));
    double B = std::sqrt(sqr(A) + sqr(2.0 * sqr(n2) * k));
    double U = std::sqrt((A + B) / 2.0);
    double V = std::sqrt((B - A) / 2.0);

    FresnelResult result;
    result.Rs = (sqr(n1 * ct1 - U) + sqr(V)) / (sqr(n1 * ct1 + U) + sqr(V));
    result.phiS = std::atan2(2.0 * n1 * V * ct1, sqr(U) + sqr(V) - sqr(n1 * ct1)) + PI_D;

    result.Rp =
        (sqr(sqr(n2) * (1.0 - sqr(k)) * ct1 - n1 * U) + sqr(2.0 * sqr(n2) * k * ct1 - n1 * V)) /
        (sqr(sqr(n2) * (1.0 - sqr(k)) * ct1 + n1 * U) + sqr(2.0 * sqr(n2) * k * ct1 + n1 * V));
    result.phiP = std::atan2(2.0 * n1 * sqr(n2) * ct1 * (2.0 * k * U - (1.0 - sqr(k)) * V),
                             sqr(sqr(n2) * (1.0 + sqr(k)) * ct1) - sqr(n1) * (sqr(U) + sqr(V)));

    return result;
}

/*
 * Evaluates the Fourier transform of the CIE XYZ sensitivity curves using Gaussian fits,
 * with "opd" (optical path difference) in micrometers and "shift" in radians.
 */
Vec3 evalSensitivity(double opd, double shift)
{
    const Vec3 val(5.4856e-13, 4.4201e-13, 5.2481e-13);
    const Vec3 pos(1.6810e+06, 1.7953e+06, 2.2084e+06);
    const Vec3 var(4.3278e+09, 9.3046e+09, 6.6121e+09);

    double phase = 2.0 * PI_D * opd * 1.0e-6;

    Vec3 xyz;
    for (int i = 0; i < 3; ++i) {
        xyz[i] = val[i] * std::sqrt(2.0 * PI_D * var[i]) * std::cos(pos[i] * phase + shift) *
                 std::exp(-var[i] * phase * phase);
    }

    // An additional Gaussian term that improves the fit of the X curve.
    xyz[0] += 9.7470e-14 * std::sqrt(2.0 * PI_D * 4.5282e+09) *
              std::cos(2.2399e+06 * phase + shift) * std::exp(-4.5282e+09 * phase * phase);

    return xyz / 1.0685e-7;
}

} // namespace

Vec3 VaryingIridescenceMicrofacet::compute(const Vec3& L,
                                           const Vec3& V,
                                           const Vec3& N,
                                           double      roughness,
                                           double      filmThickness,
                                           double      filmIor,
                                           double      baseIor,
                                           double      baseExtinctionCoefficient)
{
    double dotNL = N.dot(L);
    double dotNV = N.dot(V);
    if (dotNL <= 0.0 || dotNV <= 0.0) {
        return Vec3::Zero();
    }

    Vec3   H = (L + V).normalized();
    double dotNH = N.dot(H);
    double cosTheta1 = clamp(H.dot(L), -1.0, 1.0);

    // Convert from nanometers to micrometers.
    double thickness = filmThickness * 1.0e-3;

    // The optical path difference at normal incidence.
    double Dinc = 2.0 * filmIor * thickness;

    // Force the film's IOR to 1 (i.e., no film) as the thickness approaches zero to avoid a
    // discontinuity.
    double filmIorEff = lerp(1.0, filmIor, smoothstep(0.0, 0.03, Dinc));

    double sinTheta1Sq = 1.0 - sqr(cosTheta1);
    double cosTheta2 = std::sqrt(std::max(1.0 - sinTheta1Sq / sqr(filmIorEff), 0.0));

    // Reflectance and phase shift at the external medium/film interface.
    FresnelResult f12 = fresnelDielectric(cosTheta1, 1.0, filmIorEff);
    double        T121s = 1.0 - f12.Rs;
    double        T121p = 1.0 - f12.Rp;
    double        phi21s = PI_D - f12.phiS;
    double        phi21p = PI_D - f12.phiP;

    // Reflectance and phase shift at the film/base interface.
    FresnelResult f23 = fresnelConductor(cosTheta2, filmIorEff, baseIor, baseExtinctionCoefficient);

    // Optical path difference and the total phase shift of one round trip in the film.
    double opd = Dinc * cosTheta2;
    double phi2s = phi21s + f23.phiS;
    double phi2p = phi21p + f23.phiP;

    double r123s = std::sqrt(f12.Rs * f23.Rs);
    double r123p = std::sqrt(f12.Rp * f23.Rp);

    // Reflectance term for m = 0 (the DC term amplitude).
    double Rs = sqr(T121s) * f23.Rs / (1.0 - f12.Rs * f23.Rs);
    double Rp = sqr(T121p) * f23.Rp / (1.0 - f12.Rp * f23.Rp);
    Vec3   xyz = 0.5 * (f12.Rs + Rs + f12.Rp + Rp) * evalSensitivity(0.0, 0.0);

    // Reflectance terms for m > 0 (pairs of Dirac deltas), truncated to three terms.
    double Cms = Rs - T121s;
    double Cmp = Rp - T121p;
    for (int m = 1; m <= 3; ++m) {
        Cms *= r123s;
        Cmp *= r123p;
        Vec3 Sms = 2.0 * evalSensitivity(m * opd, m * phi2s);
        Vec3 Smp = 2.0 * evalSensitivity(m * opd, m * phi2p);
        xyz += 0.5 * (Cms * Sms + Cmp * Smp);
    }

    // CIE XYZ (with the E illuminant as the white point) to linear RGB.
    static const Eigen::Matrix3d xyzToRgb =
        (Eigen::Matrix3d() <<
            2.3706743, -0.9000405, -0.4706338,
            -0.5138850, 1.4253036, 0.0885814,
            0.0052982, -0.0146949, 1.0093968).finished();
    Vec3 specularColor = (xyzToRgb * xyz).cwiseMax(0.0).cwiseMin(1.0);

    double alpha = roughness * roughness;
    double sqAlpha = alpha * alpha;
    double G = Ggx::computeG1(dotNL, sqAlpha) * Ggx::computeG1(dotNV, sqAlpha);
    double D = Ggx::computeD(dotNH, sqAlpha);

    return specularColor * D * G / (4.0 * dotNL * dotNV);
}
