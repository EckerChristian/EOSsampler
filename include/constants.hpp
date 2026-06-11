#pragma once

namespace eos {

inline constexpr double G = 6.674E-11;
inline constexpr double Gh = 0.17229548681785437;
inline constexpr double c = 299792458.0;

inline constexpr double M_solar = 1.988E30;
inline constexpr double M_solarh = 8.572849497150262;

inline constexpr double pi = 3.14159265358979323846;

inline constexpr double mtoSI = (c * c) / G;
inline constexpr double ptoSI = (c * c * c * c) / G;

inline constexpr double gevtoinvfm = 5.067730756672362;
inline constexpr double MeVbyfm3ToGev4 = 7.683505392094774E-6;

inline constexpr double ns = 0.16;

// pQCD constants
inline constexpr double a11 = -2.0;
inline constexpr double a21 = -1.0;
inline constexpr double a22 = -3.0;
inline constexpr double a23 = -5.0021;
inline constexpr double a3 = 1.67026;

inline constexpr double nf = 3.0;
inline constexpr double nc = 3.0;
inline constexpr double da = 8.0;
inline constexpr double ca = nc;
inline constexpr double cf = da / 2.0 / nc;
inline constexpr double tf = 0.5;

inline constexpr double beta0 = 11.0 / 3.0 * ca - 4.0 / 3.0 * tf * nf;
inline constexpr double beta1 =
    34.0 / 3.0 * ca * ca - 4.0 * cf * tf * nf - 20.0 / 3.0 * ca * tf * nf;

inline constexpr double beta2 =
    2857.0 / 54.0 * ca * ca * ca
    + 2.0 * cf * cf * tf * nf
    - 205.0 / 9.0 * cf * ca * tf * nf
    - 1415.0 / 27.0 * ca * ca * tf * nf
    + 44.0 / 9.0 * cf * tf * tf * nf * nf
    + 158.0 / 27.0 * ca * tf * tf * nf * nf;

inline constexpr double LMSbar = 0.34408323550919756;

inline constexpr double Xmin = 0.5;
inline constexpr double Xmax = 2.0;

} // namespace eos
