/*
 * The implementation of "Multiple-Scattering Microfacet BSDFs with the Smith Model".
 * The original code is available at:
 * https://eheitzresearch.wordpress.com/240-2/
 *
 * OpenGL Mathematics (GLM) functions in the original code are replaced by Eigen.
 */

#include <libbsdf/ReflectanceModel/MultipleScatteringSmith.h>

using namespace lb;

/* Abstract class for height distribution */
class MicrosurfaceHeight
{
public:
    virtual ~MicrosurfaceHeight() {}

    // height PDF
    virtual double P1(double h) const = 0;

    // height CDF
    virtual double C1(double h) const = 0;

    // inverse of the height CDF
    virtual double invC1(double U) const = 0;
};

/* Uniform height distribution in [-1, 1] */
class MicrosurfaceHeightUniform : public MicrosurfaceHeight
{
public:
    virtual ~MicrosurfaceHeightUniform() {}

    // height PDF
    virtual double P1(double h) const;

    // height CDF
    virtual double C1(double h) const;

    // inverse of the height CDF
    virtual double invC1(double U) const;
};

/* Gaussian height distribution N(0,1) */
class MicrosurfaceHeightGaussian : public MicrosurfaceHeight
{
public:
    virtual ~MicrosurfaceHeightGaussian() {}

    // height PDF
    virtual double P1(double h) const;

    // height CDF
    virtual double C1(double h) const;

    // inverse of the height CDF
    virtual double invC1(double U) const;
};

/* Abstract class for slope distribution */
class MicrosurfaceSlope
{
public:
    MicrosurfaceSlope(double alpha_x = 1, double alpha_y = 1)
        : m_alpha_x(alpha_x), m_alpha_y(alpha_y)
    {
    }

    virtual ~MicrosurfaceSlope() {}

    // roughness
    double m_alpha_x, m_alpha_y;

    // projected roughness in wi
    double alpha_i(const Vec3& wi) const;

    // distribution of normals (NDF)
    double D(const Vec3& wm) const;

    // distribution of visible normals (VNDF)
    double D_wi(const Vec3& wi, const Vec3& wm) const;

    // sample the VNDF
    Vec3 sampleD_wi(const Vec3& wi, double U1, double U2) const;

    // distribution of slopes
    virtual double P22(double slope_x, double slope_y) const = 0;

    // Smith's Lambda function
    virtual double Lambda(const Vec3& wi) const = 0;

    // projected area towards incident direction
    virtual double projectedArea(const Vec3& wi) const = 0;

    // sample the distribution of visible slopes with alpha=1.0
    virtual Vec2 sampleP22_11(double theta_i, double U1, double U2) const = 0;
};

/* Beckmann slope distribution */
class MicrosurfaceSlopeBeckmann : public MicrosurfaceSlope
{
public:
    MicrosurfaceSlopeBeckmann(double alpha_x = 1, double alpha_y = 1)
        : MicrosurfaceSlope(alpha_x, alpha_y)
    {
    }

    virtual ~MicrosurfaceSlopeBeckmann() {}

    // distribution of slopes
    virtual double P22(double slope_x, double slope_y) const;

    // Smith's Lambda function
    virtual double Lambda(const Vec3& wi) const;

    // projected area towards incident direction
    virtual double projectedArea(const Vec3& wi) const;

    // sample the distribution of visible slopes with alpha=1.0
    virtual Vec2 sampleP22_11(double theta_i, double U1, double U2) const;
};

/* GGX slope distribution */
class MicrosurfaceSlopeGGX : public MicrosurfaceSlope
{
public:
    MicrosurfaceSlopeGGX(double alpha_x = 1, double alpha_y = 1)
        : MicrosurfaceSlope(alpha_x, alpha_y) {}

    virtual ~MicrosurfaceSlopeGGX() {}

    // distribution of slopes
    virtual double P22(double slope_x, double slope_y) const;

    // Smith's Lambda function
    virtual double Lambda(const Vec3& wi) const;

    // projected area towards incident direction
    virtual double projectedArea(const Vec3& wi) const;

    // sample the distribution of visible slopes with alpha=1.0
    virtual Vec2 sampleP22_11(double theta_i, double U1, double U2) const;
};

/* Abstract class for microsurface */
class Microsurface
{
public:
    // height distribution
    const MicrosurfaceHeight* m_microsurfaceheight;

    // slope distribution
    const MicrosurfaceSlope* m_microsurfaceslope;

    Microsurface(bool   height_uniform, // uniform or Gaussian height distribution
                 bool   slope_beckmann, // Beckmann or GGX slope distribution
                 double alpha_x,
                 double alpha_y)
                 : m_microsurfaceheight((height_uniform) ? static_cast<MicrosurfaceHeight*>(new MicrosurfaceHeightUniform)
                                                         : static_cast<MicrosurfaceHeight*>(new MicrosurfaceHeightGaussian)),
                   m_microsurfaceslope((slope_beckmann) ? static_cast<MicrosurfaceSlope*>(new MicrosurfaceSlopeBeckmann(alpha_x, alpha_y))
                                                        : static_cast<MicrosurfaceSlope*>(new MicrosurfaceSlopeGGX(alpha_x, alpha_y))) {}

    virtual ~Microsurface()
    {
        delete m_microsurfaceheight;
        delete m_microsurfaceslope;
    }

    // evaluate BSDF with a random walk (stochastic but unbiased)
    // scatteringOrder=0 --> contribution from all scattering events
    // scatteringOrder=1 --> contribution from 1st bounce only
    // scatteringOrder=2 --> contribution from 2nd bounce only, etc..
    virtual double eval(const Vec3& wi, const Vec3& wo, int scatteringOrder = 0) const;

    // sample BSDF with a random walk
    // scatteringOrder is set to the number of bounces computed for this sample
    virtual Vec3 sample(const Vec3& wi, int& scatteringOrder) const;
    Vec3 sample(const Vec3& wi) const
    {
        int scatteringOrder;
        return sample(wi, scatteringOrder);
    }

    // masking function
    double G_1(const Vec3& wi) const;

    // masking function at height h0
    double G_1(const Vec3& wi, double h0) const;

    // sample height in outgoing direction
    double sampleHeight(const Vec3& wo, double h0, double U) const;

    // evaluate local phase function
    virtual double evalPhaseFunction(const Vec3& wi, const Vec3& wo) const = 0;

    // sample local phase function
    virtual Vec3 samplePhaseFunction(const Vec3& wi) const = 0;

    // evaluate BSDF limited to single scattering
    // this is in average equivalent to eval(wi, wo, 1);
    virtual double evalSingleScattering(const Vec3& wi, const Vec3& wo) const = 0;
};

/* Microsurface made of conductor material */
class MicrosurfaceConductor : public Microsurface
{
public:
    MicrosurfaceConductor(bool   height_uniform, // uniform or Gaussian
                          bool   slope_beckmann, // Beckmann or GGX
                          double alpha_x,
                          double alpha_y)
        : Microsurface(height_uniform, slope_beckmann, alpha_x, alpha_y)
    {
    }

    virtual ~MicrosurfaceConductor() {}

    // evaluate local phase function
    virtual double evalPhaseFunction(const Vec3& wi, const Vec3& wo) const;

    // sample local phase function
    virtual Vec3 samplePhaseFunction(const Vec3& wi) const;

    // evaluate BSDF limited to single scattering
    // this is in average equivalent to eval(wi, wo, 1);
    virtual double evalSingleScattering(const Vec3& wi, const Vec3& wo) const;
};

/* Microsurface made of conductor material */
class MicrosurfaceDielectric : public Microsurface
{
public:
    MicrosurfaceDielectric(bool   height_uniform, // uniform or Gaussian
                           bool   slope_beckmann, // Beckmann or GGX
                           double alpha_x,
                           double alpha_y,
                           double eta = 1.5)
        : Microsurface(height_uniform, slope_beckmann, alpha_x, alpha_y), m_eta(eta)
    {
    }

    virtual ~MicrosurfaceDielectric() {}

    // evaluate BSDF with a random walk (stochastic but unbiased)
    // scatteringOrder=0 --> contribution from all scattering events
    // scatteringOrder=1 --> contribution from 1st bounce only
    // scatteringOrder=2 --> contribution from 2nd bounce only, etc..
    virtual double eval(const Vec3& wi, const Vec3& wo, int scatteringOrder = 0) const;

    // sample final BSDF with a random walk
    // scatteringOrder is set to the number of bounces computed for this sample
    virtual Vec3 sample(const Vec3& wi, int& scatteringOrder) const;

    // evaluate local phase function
    virtual double evalPhaseFunction(const Vec3& wi, const Vec3& wo) const;
    double
    evalPhaseFunction(const Vec3& wi, const Vec3& wo, bool wi_outside, bool wo_outside) const;

    // sample local phase function
    virtual Vec3 samplePhaseFunction(const Vec3& wi) const;
    Vec3         samplePhaseFunction(const Vec3& wi, bool wi_outside, bool& wo_outside) const;

    // evaluate BSDF limited to single scattering
    // this is in average equivalent to eval(wi, wo, 1);
    virtual double evalSingleScattering(const Vec3& wi, const Vec3& wo) const;

private:
    double fresnel(const Vec3& wi, const Vec3& wm, double eta) const;
    Vec3   refract(const Vec3& wi, const Vec3& wm, double eta) const;

    double m_eta;
};

/* Microsurface made of conductor material */
class MicrosurfaceDiffuse : public Microsurface
{
public:
    MicrosurfaceDiffuse(bool   height_uniform, // uniform or Gaussian
                        bool   slope_beckmann, // Beckmann or GGX
                        double alpha_x,
                        double alpha_y)
                        : Microsurface(height_uniform, slope_beckmann, alpha_x, alpha_y) {}

    virtual ~MicrosurfaceDiffuse() {}

    // evaluate local phase function
    virtual double evalPhaseFunction(const Vec3& wi, const Vec3& wo) const;

    // sample local phase function
    virtual Vec3 samplePhaseFunction(const Vec3& wi) const;

    // evaluate BSDF limited to single scattering
    // this is in average equivalent to eval(wi, wo, 1);
    virtual double evalSingleScattering(const Vec3& wi, const Vec3& wo) const;
};

//
// Implementation
//

#include <cfloat>
#include <cmath>
#include <algorithm>
#include <random>
using namespace std;

#ifndef M_PI
#define M_PI            3.14159265358979323846 /* pi */
#endif

#ifndef M_PI_2
#define M_PI_2          1.57079632679489661923 /* pi/2 */
#endif

#define INV_M_PI        0.31830988618379067153 /* 1/pi */
#define SQRT_M_PI       1.77245385090551602729 /* sqrt(pi) */
#define SQRT_2          1.41421356237309504880 /* sqrt(2) */
#define INV_SQRT_M_PI   0.56418958354775628694 /* 1/sqrt(pi) */
#define INV_2_SQRT_M_PI 0.28209479177387814347 /* 0.5/sqrt(pi) */
#define INV_SQRT_2_M_PI 0.3989422804014326779  /* 1/sqrt(2*pi) */
#define INV_SQRT_2      0.7071067811865475244  /* 1/sqrt(2) */

#include <libbsdf/Common/Utility.h>

static bool IsFiniteNumber(double x) { return (x <= FLT_MAX && x >= -FLT_MAX); }

static std::random_device rd;
static std::mt19937 gen(rd());
static std::uniform_real_distribution<> dis(0, 1);

static double generateRandomNumber() { return dis(gen); }

static double erfAs(double x)
{
    // constants
    double a1 = 0.254829592;
    double a2 = -0.284496736;
    double a3 = 1.421413741;
    double a4 = -1.453152027;
    double a5 = 1.061405429;
    double p = 0.3275911;

    // Save the sign of x
    int sign = 1;
    if (x < 0) sign = -1;
    x = fabs(x);

    // A&S formula 7.1.26
    double t = 1.0 / (1.0 + p * x);
    double y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * exp(-x * x);

    return sign*y;
}

static double erfinv(double x)
{
    double w, p;
    w = -log((1 - x) * (1 + x));

    if (w < 5) {
        w = w - 2.5;
        p = 2.81022636e-08;
        p = 3.43273939e-07 + p * w;
        p = -3.5233877e-06 + p * w;
        p = -4.39150654e-06 + p * w;
        p = 0.00021858087 + p * w;
        p = -0.00125372503 + p * w;
        p = -0.00417768164 + p * w;
        p = 0.246640727 + p * w;
        p = 1.50140941 + p * w;
    }
    else {
        w = sqrt(w) - 3;
        p = -0.000200214257;
        p = 0.000100950558 + p * w;
        p = 0.00134934322 + p * w;
        p = -0.00367342844 + p * w;
        p = 0.00573950773 + p * w;
        p = -0.0076224613 + p * w;
        p = 0.00943887047 + p * w;
        p = 1.00167406 + p * w;
        p = 2.83297682 + p * w;
    }

    return p * x;
}

/* A method to compute the gamma() function. */
static double abgam(double x)
{
    double gam[10];

    gam[0] = 1. / 12.;
    gam[1] = 1. / 30.;
    gam[2] = 53. / 210.;
    gam[3] = 195. / 371.;
    gam[4] = 22999. / 22737.;
    gam[5] = 29944523. / 19733142.;
    gam[6] = 109535241009. / 48264275462.;

    return  0.5 * log(2 * M_PI) - x + (x - 0.5) * log(x)
            + gam[0] / (x + gam[1] / (x + gam[2] / (x + gam[3] / (x + gam[4] / (x + gam[5] / (x + gam[6] / x))))));
}

namespace mss {
static double gamma(double x)
{
    return exp(abgam(x + 5)) / (x * (x + 1) * (x + 2) * (x + 3) * (x + 4));
}

static double beta(double m, double n)
{
    return (mss::gamma(m) * mss::gamma(n) / mss::gamma(m + n));
}
} // namespace mss

/*
 * Implementation of MicrosurfaceHeightUniform
 */

double MicrosurfaceHeightUniform::P1(double h) const
{
    return (h >= -1 && h <= 1) ? 0.5 : 0;
}

double MicrosurfaceHeightUniform::C1(double h) const
{
    using std::max;
    using std::min;

    return min(1.0, max(0.0, 0.5 * (h + 1)));
}

double MicrosurfaceHeightUniform::invC1(double U) const
{
    using std::max;
    using std::min;

    return max(-1.0, min(1.0, 2 * U - 1));
}

double MicrosurfaceHeightGaussian::P1(double h) const
{
    return INV_SQRT_2_M_PI * exp(-0.5 * h * h);
}

double MicrosurfaceHeightGaussian::C1(double h) const
{
    return 0.5 + 0.5 * erfAs(INV_SQRT_2 * h);
}

double MicrosurfaceHeightGaussian::invC1(double U) const
{
    return SQRT_2 * erfinv(2 * U - 1);
}

/*
 * Implementation of MicrosurfaceSlope
 */

double MicrosurfaceSlope::D(const Vec3& wm) const
{
    if (wm.z() <= 0)
        return 0;

    // slope of wm
    double slope_x = -wm.x() / wm.z();
    double slope_y = -wm.y() / wm.z();

    return P22(slope_x, slope_y) / (wm.z() * wm.z() * wm.z() * wm.z());
}

double MicrosurfaceSlope::D_wi(const Vec3& wi, const Vec3& wm) const
{
    if (wm.z() <= 0)
        return 0;

    // normalization coefficient
    double projectedarea = projectedArea(wi);
    if (projectedarea == 0)
        return 0;

    double c = 1 / projectedarea;

    using std::max;

    return c * max(0.0, wi.dot(wm)) * D(wm);
}

Vec3 MicrosurfaceSlope::sampleD_wi(const Vec3& wi, double U1, double U2) const
{
    // stretch to match configuration with alpha=1.0
    Vec3 wi_11 = Vec3(m_alpha_x * wi.x(), m_alpha_y * wi.y(), wi.z()).normalized();

    // sample visible slope with alpha=1.0
    Vec2 slope_11 = sampleP22_11(acos(wi_11.z()), U1, U2);

    // align with view direction
    double phi = atan2(wi_11.y(), wi_11.x());
    Vec2 slope(cos(phi) * slope_11.x() - sin(phi) * slope_11.y(),
                sin(phi) * slope_11.x() + cos(phi) * slope_11.y());

    // stretch back
    slope.x() *= m_alpha_x;
    slope.y() *= m_alpha_y;

    // if numerical instability
    if (slope.x() != slope.x() || !IsFiniteNumber(slope.x())) {
        if (wi.z() > 0) {
            return Vec3(0, 0, 1);
        }
        else {
            return Vec3(wi.x(), wi.y(), 0).normalized();
        }
    }

    // compute normal
    Vec3 wm = Vec3(-slope.x(), -slope.y(), 1).normalized();
    return wm;
}

double MicrosurfaceSlope::alpha_i(const Vec3& wi) const
{
    double invSinTheta2 = 1 / (1 - wi.z() * wi.z());
    double cosPhi2 = wi.x() * wi.x() * invSinTheta2;
    double sinPhi2 = wi.y() * wi.y() * invSinTheta2;
    double alpha_i = sqrt(cosPhi2 * m_alpha_x * m_alpha_x + sinPhi2 * m_alpha_y * m_alpha_y);
    return alpha_i;
}

double MicrosurfaceSlopeBeckmann::P22(double slope_x, double slope_y) const
{
    return 1 / (M_PI * m_alpha_x * m_alpha_y) *
           exp(-slope_x * slope_x / (m_alpha_x * m_alpha_x) -
                slope_y * slope_y / (m_alpha_y * m_alpha_y));
}

double MicrosurfaceSlopeBeckmann::Lambda(const Vec3& wi) const
{
    if (wi.z() > 0.9999)
        return 0;

    if (wi.z() < -0.9999)
        return -1;

    // a
    double theta_i = acos(wi.z());
    double a = 1 / tan(theta_i) / alpha_i(wi);

    return 0.5 * (erfAs(a) - 1.0) + INV_2_SQRT_M_PI / a * exp(-a * a);
}

double MicrosurfaceSlopeBeckmann::projectedArea(const Vec3& wi) const
{
    if (wi.z() > 0.9999)
        return 1;

    if (wi.z() < -0.9999)
        return 0;

    // a
    double alphai = alpha_i(wi);
    double theta_i = acos(wi.z());
    double a = 1 / tan(theta_i) / alphai;

    return 0.5 * (erfAs(a) + 1.0) * wi.z() + INV_2_SQRT_M_PI * alphai * sin(theta_i) * exp(-a * a);
}

Vec2 MicrosurfaceSlopeBeckmann::sampleP22_11(double theta_i, double U, double U_2) const
{
    using std::max;
    using std::min;

    Vec2 slope2D;

    if (theta_i < 0.0001) {
        double r = sqrt(-log(U));
        double phi = 6.28318530718 * U_2;
        slope2D.x() = r * cos(phi);
        slope2D.y() = r * sin(phi);
        return slope2D;
    }

    // constant
    double sin_theta_i = sin(theta_i);
    double cos_theta_i = cos(theta_i);

    // slope associated to theta_i
    double slope_i = cos_theta_i / sin_theta_i;

    // projected area
    double a = cos_theta_i / sin_theta_i;
    double projectedarea =
        0.5 * (erfAs(a) + 1) * cos_theta_i + INV_2_SQRT_M_PI * sin_theta_i * exp(-a * a);
    if (projectedarea < 0.0001 || projectedarea != projectedarea) {
        return Vec2(0, 0);
    }

    // VNDF normalization factor
    double c = 1 / projectedarea;

    // search
    double erf_min = -0.9999;
    double erf_max = max(erf_min, erfAs(slope_i));
    double erf_current = 0.5 * (erf_min + erf_max);

    while (erf_max - erf_min > 0.00001) {
        if (!(erf_current >= erf_min && erf_current <= erf_max)) {
            erf_current = 0.5 * (erf_min + erf_max);
        }

        // evaluate slope
        double slope = erfinv(erf_current);

        // CDF
        double CDF = (slope >= slope_i) ? 1
                                        : c * (INV_2_SQRT_M_PI * sin_theta_i * exp(-slope * slope) +
                                               cos_theta_i * (0.5 + 0.5 * erfAs(slope)));
        double diff = CDF - U;

        // test estimate
        if (abs(diff) < 0.00001)
            break;

        // update bounds
        if (diff > 0) {
            if (erf_max == erf_current)
                break;

            erf_max = erf_current;
        }
        else {
            if (erf_min == erf_current)
                break;

            erf_min = erf_current;
        }

        // update estimate
        double derivative = 0.5 * c * cos_theta_i - 0.5 * c * sin_theta_i * slope;
        erf_current -= diff / derivative;
    }

    slope2D.x() = erfinv(min(erf_max, max(erf_min, erf_current)));
    slope2D.y() = erfinv(2 * U_2 - 1);

    return slope2D;
}

double MicrosurfaceSlopeGGX::P22(double slope_x, double slope_y) const
{
    double temp = 1 + slope_x * slope_x / (m_alpha_x * m_alpha_x) +
                  slope_y * slope_y / (m_alpha_y * m_alpha_y);

    return 1 / (M_PI * m_alpha_x * m_alpha_y) / (temp * temp);
}

double MicrosurfaceSlopeGGX::Lambda(const Vec3& wi) const
{
    if (wi.z() > 0.9999)
        return 0;

    if (wi.z() < -0.9999)
        return -1;

    // a
    double theta_i = acos(wi.z());
    double a = 1 / tan(theta_i) / alpha_i(wi);

    return 0.5 * (-1 + sign(a) * sqrt(1 + 1 / (a * a)));
}

double MicrosurfaceSlopeGGX::projectedArea(const Vec3& wi) const
{
    if (wi.z() > 0.9999)
        return 1;

    if (wi.z() < -0.9999)
        return 0;

    // a
    double theta_i = acos(wi.z());
    double sin_theta_i = sin(theta_i);

    double alphai = alpha_i(wi);

    return 0.5 * (wi.z() + sqrt(wi.z() * wi.z() + sin_theta_i * sin_theta_i * alphai * alphai));
}

Vec2 MicrosurfaceSlopeGGX::sampleP22_11(double theta_i, double U, double U_2) const
{
    Vec2 slope;

    if (theta_i < 0.0001) {
        double r = sqrt(U / (1 - U));
        double phi = 6.28318530718 * U_2;
        slope.x() = r * cos(phi);
        slope.y() = r * sin(phi);

        return slope;
    }

    // constant
    double sin_theta_i = sin(theta_i);
    double cos_theta_i = cos(theta_i);
    double tan_theta_i = sin_theta_i / cos_theta_i;

    // slope associated to theta_i
    //double slope_i = cos_theta_i / sin_theta_i;

    // projected area
    double projectedarea = 0.5 * (cos_theta_i + 1);
    if (projectedarea < 0.0001 || projectedarea != projectedarea) {
        return Vec2(0, 0);
    }
    // normalization coefficient
    double c = 1 / projectedarea;

    double A = 2 * U / cos_theta_i / c - 1;
    double B = tan_theta_i;
    double tmp = 1 / (A * A - 1);

    using std::max;

    double D = sqrt(max(0.0, B * B * tmp * tmp - (A * A - B * B) * tmp));
    double slope_x_1 = B * tmp - D;
    double slope_x_2 = B * tmp + D;
    slope.x() = (A < 0 || slope_x_2 > 1 / tan_theta_i) ? slope_x_1 : slope_x_2;

    double U2;
    double S;
    if (U_2 > 0.5) {
        S = 1;
        U2 = 2 * (U_2 - 0.5);
    }
    else {
        S = -1;
        U2 = 2 * (0.5 - U_2);
    }
    double z = (U2 * (U2 * (U2 * 0.27385 - 0.73369) + 0.46341)) /
               (U2 * (U2 * (U2 * 0.093073 + 0.309420) - 1.000000) + 0.597999);
    slope.y() = S * z * sqrt(1 + slope.x() * slope.x());

    return slope;
}

/*
 * Implementation of Microsurface
 */

double Microsurface::G_1(const Vec3& wi) const
{
    if (wi.z() > 0.9999)
        return 1;

    if (wi.z() <= 0)
        return 0;

    // Lambda
    double Lambda = m_microsurfaceslope->Lambda(wi);

    return 1 / (1 + Lambda);
}

double Microsurface::G_1(const Vec3& wi, double h0) const
{
    if (wi.z() > 0.9999)
        return 1;

    if (wi.z() <= 0)
        return 0;

    // height CDF
    double C1_h0 = m_microsurfaceheight->C1(h0);
    // Lambda
    double Lambda = m_microsurfaceslope->Lambda(wi);

    return pow(C1_h0, Lambda);
}

double Microsurface::sampleHeight(const Vec3& wr, double hr, double U) const
{
    if (wr.z() > 0.9999)
        return FLT_MAX;

    if (wr.z() < -0.9999) {
        return m_microsurfaceheight->invC1(U * m_microsurfaceheight->C1(hr));
    }

    if (fabs(wr.z()) < 0.0001)
        return hr;

    // probability of intersection
    double G_1_ = G_1(wr, hr);

    if (U > 1 - G_1_) {
        // leave the microsurface
        return FLT_MAX;
    }

    double h = m_microsurfaceheight->invC1(
        m_microsurfaceheight->C1(hr) / pow((1 - U), 1 / m_microsurfaceslope->Lambda(wr)));
    return h;
}

Vec3 Microsurface::sample(const Vec3& wi, int& scatteringOrder) const
{
    // init
    Vec3 wr = -wi;
    double hr = 1 + m_microsurfaceheight->invC1(0.999);

    // random walk
    scatteringOrder = 0;
    while (true) {
        // next height
        double U = generateRandomNumber();
        hr = sampleHeight(wr, hr, U);

        // leave the microsurface?
        if (hr == FLT_MAX) {
            break;
        }
        else {
            scatteringOrder++;
        }

        // next direction
        wr = samplePhaseFunction(-wr);

        // if NaN (should not happen, just in case)
        if ((hr != hr) || (wr.z() != wr.z())) {
            return Vec3(0, 0, 1);
        }
    }

    return wr;
}

double Microsurface::eval(const Vec3& wi, const Vec3& wo, int scatteringOrder) const
{
    if (wo.z() < 0)
        return 0;

    // init
    Vec3 wr = -wi;
    double hr = 1 + m_microsurfaceheight->invC1(0.999);

    double sum = 0;

    // random walk
    int current_scatteringOrder = 0;
    while (scatteringOrder == 0 || current_scatteringOrder <= scatteringOrder) {
        // next height
        double U = generateRandomNumber();
        hr = sampleHeight(wr, hr, U);

        // leave the microsurface?
        if (hr == FLT_MAX) {
            break;
        }
        else {
            current_scatteringOrder++;
        }

        // next event estimation
        double phasefunction = evalPhaseFunction(-wr, wo);
        double shadowing = G_1(wo, hr);
        double I = phasefunction * shadowing;

        if (IsFiniteNumber(I) &&
            (scatteringOrder == 0 || current_scatteringOrder == scatteringOrder)) {
            sum += I;
        }

        // next direction
        wr = samplePhaseFunction(-wr);

        // if NaN (should not happen, just in case)
        if ((hr != hr) || (wr.z() != wr.z())) {
            return 0;
        }
    }

    return sum;
}

double MicrosurfaceConductor::evalPhaseFunction(const Vec3& wi, const Vec3& wo) const
{
    // half vector
    Vec3 wh = (wi + wo).normalized();
    if (wh.z() < 0)
        return 0;

    return 0.25 * m_microsurfaceslope->D_wi(wi, wh) / wi.dot(wh);
}

Vec3 MicrosurfaceConductor::samplePhaseFunction(const Vec3& wi) const
{
    double U1 = generateRandomNumber();
    double U2 = generateRandomNumber();

    Vec3 wm = m_microsurfaceslope->sampleD_wi(wi, U1, U2);

    // reflect
    Vec3 wo = -wi + 2 * wm * wi.dot(wm);

    return wo;
}

double MicrosurfaceConductor::evalSingleScattering(const Vec3& wi, const Vec3& wo) const
{
    // half-vector
    Vec3 wh = (wi + wo).normalized();
    double D = m_microsurfaceslope->D(wh);

    // masking-shadowing
    double G2 = 1 / (1 + m_microsurfaceslope->Lambda(wi) + m_microsurfaceslope->Lambda(wo));

    // BRDF * cos
    return D * G2 / (4 * wi.z());
}

Vec3 MicrosurfaceDielectric::refract(const Vec3& wi, const Vec3& wm, double eta) const
{
    using std::max;

    double cos_theta_i = wi.dot(wm);
    double cos_theta_t2 = 1 - (1 - cos_theta_i * cos_theta_i) / (eta * eta);
    double cos_theta_t = -sqrt(max(0.0, cos_theta_t2));

    return wm * (wi.dot(wm) / eta + cos_theta_t) - wi / eta;
}

double MicrosurfaceDielectric::fresnel(const Vec3& wi, const Vec3& wm, double eta) const
{
    double cos_theta_i = wi.dot(wm);
    double cos_theta_t2 = 1 - (1 - cos_theta_i * cos_theta_i) / (eta * eta);

    // total internal reflection
    if (cos_theta_t2 <= 0)
        return 1;

    double cos_theta_t = sqrt(cos_theta_t2);

    double Rs = (cos_theta_i - eta * cos_theta_t) / (cos_theta_i + eta * cos_theta_t);
    double Rp = (eta * cos_theta_i - cos_theta_t) / (eta * cos_theta_i + cos_theta_t);

    double F = 0.5 * (Rs * Rs + Rp * Rp);
    return F;
}

// wrapper (only for the API and testing)
double MicrosurfaceDielectric::evalPhaseFunction(const Vec3& wi, const Vec3& wo) const
{
    return evalPhaseFunction(wi, wo, true, true) +
           evalPhaseFunction(wi, wo, true, false);
}

double MicrosurfaceDielectric::evalPhaseFunction(const Vec3& wi,
                                                 const Vec3& wo,
                                                 bool        wi_outside,
                                                 bool        wo_outside) const
{
    double eta = wi_outside ? m_eta : 1 / m_eta;

    if (wi_outside == wo_outside) { // reflection
        // half vector
        Vec3 wh = (wi + wo).normalized();

        return wi_outside
                 ? (0.25 * m_microsurfaceslope->D_wi(wi, wh) / wi.dot(wh) * fresnel(wi, wh, eta))
                 : (0.25 * m_microsurfaceslope->D_wi(-wi, -wh) / (-wi).dot(-wh) *
                    fresnel(-wi, -wh, eta));
    }
    else { // transmission
        Vec3 wh = -((wi + wo * eta).normalized());
        wh *= wi_outside ? sign(wh.z()) : -sign(wh.z());

        if (wh.dot(wi) < 0)
            return 0;

        using std::max;

        if (wi_outside) {
            return eta * eta * (1 - fresnel(wi, wh, eta)) *
                   m_microsurfaceslope->D_wi(wi, wh) * max(0.0, -wo.dot(wh)) *
                   1 / pow(wi.dot(wh) + eta * wo.dot(wh), 2);
        }
        else {
            return eta * eta * (1 - fresnel(-wi, -wh, eta)) *
                   m_microsurfaceslope->D_wi(-wi, -wh) * max(0.0, -(-wo).dot(-wh)) *
                   1 / pow((-wi).dot(-wh) + eta * (-wo).dot(-wh), 2);
        }
    }
}

Vec3 MicrosurfaceDielectric::samplePhaseFunction(const Vec3& wi) const
{
    bool wo_outside;
    return samplePhaseFunction(wi, true, wo_outside);
}

Vec3 MicrosurfaceDielectric::samplePhaseFunction(const Vec3& wi,
                                                 bool        wi_outside,
                                                 bool&       wo_outside) const
{
    double U1 = generateRandomNumber();
    double U2 = generateRandomNumber();

    double eta = wi_outside ? m_eta : 1 / m_eta;

    Vec3 wm = wi_outside ? (m_microsurfaceslope->sampleD_wi(wi, U1, U2))
                         : (-m_microsurfaceslope->sampleD_wi(-wi, U1, U2));

    double F = fresnel(wi, wm, eta);

    if (generateRandomNumber() < F) {
        Vec3 wo = -wi + 2 * wm * wi.dot(wm); // reflect
        return wo;
    }
    else {
        wo_outside = !wi_outside;
        Vec3 wo = refract(wi, wm, eta);
        return wo.normalized();
    }
}

double MicrosurfaceDielectric::evalSingleScattering(const Vec3& wi, const Vec3& wo) const
{
    //bool wi_outside = true;
    bool wo_outside = wo.z() > 0;

    double eta = m_eta;

    if (wo_outside) { // reflection
        // D
        Vec3   wh = (wi + wo).normalized();
        double D = m_microsurfaceslope->D(wh);

        // masking shadowing
        double Lambda_i = m_microsurfaceslope->Lambda(wi);
        double Lambda_o = m_microsurfaceslope->Lambda(wo);
        double G2 = 1 / (1 + Lambda_i + Lambda_o);

        // BRDF
        return fresnel(wi, wh, eta) * D * G2 / (4 * wi.z());
    }
    else { // refraction
        // D
        Vec3 wh = -((wi + wo * eta).normalized());
        if (eta < 1) wh = -wh;
        double D = m_microsurfaceslope->D(wh);

        // G2
        double Lambda_i = m_microsurfaceslope->Lambda(wi);
        double Lambda_o = m_microsurfaceslope->Lambda(-wo);
        double G2 = mss::beta(1 + Lambda_i, 1 + Lambda_o);

        using std::max;

        // BSDF
        return max(0.0, wi.dot(wh)) * max(0.0, -wo.dot(wh)) * 1 / wi.z() * eta * eta *
               (1 - fresnel(wi, wh, eta)) * G2 * D / pow(wi.dot(wh) + eta * wo.dot(wh), 2);
    }
}

double MicrosurfaceDielectric::eval(const Vec3& wi, const Vec3& wo, int scatteringOrder) const
{
    // init
    Vec3 wr = -wi;
    double hr = 1 + m_microsurfaceheight->invC1(0.999);
    bool outside = true;

    double sum = 0;

    // random walk
    int current_scatteringOrder = 0;
    while (scatteringOrder == 0 || current_scatteringOrder <= scatteringOrder) {
        // next height
        double U = generateRandomNumber();
        hr = (outside) ? sampleHeight(wr, hr, U) : -sampleHeight(-wr, -hr, U);

        // leave the microsurface?
        if (hr == FLT_MAX || hr == -FLT_MAX) {
            break;
        }
        else {
            current_scatteringOrder++;
        }

        // next event estimation
        double phasefunction = evalPhaseFunction(-wr, wo, outside, (wo.z() > 0));
        double shadowing = (wo.z() > 0) ? G_1(wo, hr) : G_1(-wo, -hr);
        double I = phasefunction * shadowing;

        if (IsFiniteNumber(I) &&
            (scatteringOrder == 0 || current_scatteringOrder == scatteringOrder)) {
            sum += I;
        }

        // next direction
        wr = samplePhaseFunction(-wr, outside, outside);

        // if NaN (should not happen, just in case)
        if ((hr != hr) || (wr.z() != wr.z())) {
            return 0;
        }
    }

    return sum;
}

Vec3 MicrosurfaceDielectric::sample(const Vec3& wi, int& scatteringOrder) const
{
    // init
    Vec3 wr = -wi;
    double hr = 1 + m_microsurfaceheight->invC1(0.999);
    bool outside = true;

    // random walk
    scatteringOrder = 0;
    while (true) {
        // next height
        double U = generateRandomNumber();
        hr = (outside) ? sampleHeight(wr, hr, U) : -sampleHeight(-wr, -hr, U);

        // leave the microsurface?
        if (hr == FLT_MAX || hr == -FLT_MAX) {
            break;
        }
        else {
            scatteringOrder++;
        }

        // next direction
        wr = samplePhaseFunction(-wr, outside, outside);

        // if NaN (should not happen, just in case)
        if ((hr != hr) || (wr.z() != wr.z())) {
            return Vec3(0, 0, 1);
        }
    }

    return wr;
}

double MicrosurfaceDiffuse::evalPhaseFunction(const Vec3& wi, const Vec3& wo) const
{
    double U1 = generateRandomNumber();
    double U2 = generateRandomNumber();
    Vec3 wm = m_microsurfaceslope->sampleD_wi(wi, U1, U2);

    using std::max;

    return 1 / M_PI * max(0.0, wo.dot(wm));
}

// build orthonormal basis (Building an Orthonormal Basis from a 3D Unit Vector Without Normalization, [Frisvad2012])
void buildOrthonormalBasis(Vec3& omega_1, Vec3& omega_2, const Vec3& omega_3)
{
    if (omega_3.z() < -0.9999999) {
        omega_1 = Vec3(0, -1, 0);
        omega_2 = Vec3(-1, 0, 0);
    }
    else {
        double a = 1 / (1 + omega_3.z());
        double b = -omega_3.x() * omega_3.y() * a;
        omega_1 = Vec3(1 - omega_3.x() * omega_3.x() * a, b, -omega_3.x());
        omega_2 = Vec3(b, 1 - omega_3.y() * omega_3.y() * a, -omega_3.y());
    }
}

Vec3 MicrosurfaceDiffuse::samplePhaseFunction(const Vec3& wi) const
{
    double U1 = generateRandomNumber();
    double U2 = generateRandomNumber();
    double U3 = generateRandomNumber();
    double U4 = generateRandomNumber();

    Vec3 wm = m_microsurfaceslope->sampleD_wi(wi, U1, U2);

    // sample diffuse reflection
    Vec3 w1, w2;
    buildOrthonormalBasis(w1, w2, wm);

    double r1 = 2 * U3 - 1;
    double r2 = 2 * U4 - 1;

    using std::max;

    // concentric map code from
    // http://psgraphics.blogspot.ch/2011/01/improved-code-for-concentric-map.html
    double phi, r;
    if (r1 == 0 && r2 == 0) {
        r = phi = 0;
    }
    else if (r1 * r1 > r2 * r2) {
        r = r1;
        phi = (M_PI / 4) * (r2 / r1);
    }
    else {
        r = r2;
        phi = (M_PI / 2) - (r1 / r2) * (M_PI / 4);
    }
    double x = r * cos(phi);
    double y = r * sin(phi);
    double z = sqrt(max(0.0, 1 - x * x - y * y));
    Vec3 wo = x * w1 + y * w2 + z * wm;

    return wo;
}

// stochastic evaluation
// Heitz and Dupuy 2015
// Implementing a Simple Anisotropic Rough Diffuse Material with Stochastic Evaluation
double MicrosurfaceDiffuse::evalSingleScattering(const Vec3& wi, const Vec3& wo) const
{
    // sample visible microfacet
    double U1 = generateRandomNumber();
    double U2 = generateRandomNumber();
    Vec3 wm = m_microsurfaceslope->sampleD_wi(wi, U1, U2);

    // shadowing given masking
    double Lambda_i = m_microsurfaceslope->Lambda(wi);
    double Lambda_o = m_microsurfaceslope->Lambda(wo);
    double G2_given_G1 = (1 + Lambda_i) / (1 + Lambda_i + Lambda_o);

    using std::max;

    // evaluate diffuse and shadowing given masking
    return 1 / M_PI * max(0.0, wm.dot(wo)) * G2_given_G1;
}

// Interface to evaluate a BSDF value with iterations.
Vec3 MultipleScatteringSmith::compute(const Vec3&  L,
                                      const Vec3&  V,
                                      const Vec3&  color,
                                      double       roughnessX,
                                      double       roughnessY,
                                      double       refractiveIndex,
                                      MaterialType materialType,
                                      HeightType   heightType,
                                      SlopeType    slopeType,
                                      int          numIterations)
{
    bool uniformHeightUsed;
    if (heightType == UNIFORM_HEIGHT) {
        uniformHeightUsed = true;
    }
    else {
        uniformHeightUsed = false;
    }

    bool beckmannSlopeUsed;
    if (slopeType == BECKMANN_SLOPE) {
        beckmannSlopeUsed = true;
    }
    else {
        beckmannSlopeUsed = false;
    }

    double alphaX = roughnessX * roughnessX;
    double alphaY = roughnessY * roughnessY;

    Microsurface* microsurface;

    switch (materialType) {
        case CONDUCTOR_MATERIAL:
            microsurface = new MicrosurfaceConductor(
                uniformHeightUsed, beckmannSlopeUsed, alphaX, alphaY);
            break;
        case DIELECTRIC_MATERIAL:
            microsurface = new MicrosurfaceDielectric(
                uniformHeightUsed, beckmannSlopeUsed, alphaX, alphaY, refractiveIndex);
            break;
        case DIFFUSE_MATERIAL:
            microsurface = new MicrosurfaceDiffuse(
                uniformHeightUsed, beckmannSlopeUsed, alphaX, alphaY);
            break;
        default:
            lbError << "[MultipleScatteringSmith::compute] Invalid material type: " << materialType;
            return Vec3::Zero();
            break;
    }

    double sum = 0.0;
    for (int i = 0; i < numIterations; ++i) {
        Vec3 wi = L;
        Vec3 wo = V;
        sum += microsurface->eval(wi, wo, 0);
    }
    double val = sum / numIterations;

    delete microsurface;

    using std::abs;
    using std::max;

    return color * val / max(abs(V.z()), EPSILON_D);
}
