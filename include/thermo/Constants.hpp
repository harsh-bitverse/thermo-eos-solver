#ifndef THERMO_CONSTANTS_HPP
#define THERMO_CONSTANTS_HPP

#include <cmath>

namespace thermo::constants {

// ============================================================================
// Fundamental Physical & Universal Constants
// ============================================================================

/// Universal Gas Constant R in SI units [J/(mol*K)]
/// Source: CODATA 2018 / 2022 recommended value (exact defining constant).
inline constexpr double R = 8.31446261815324;

/// Standard Reference Temperature T_ref in [K] (25.0 degrees Celsius)
inline constexpr double T_REF = 298.15;

/// Standard Reference Pressure P_ref in [Pa] (1.0 bar = 100,000 Pa)
inline constexpr double P_REF = 100000.0;

/// Standard Atmospheric Pressure in [Pa] (1.0 atm = 101,325 Pa)
inline constexpr double P_ATM = 101325.0;

// ============================================================================
// Peng-Robinson (1976) Dimensionless EOS Constants
// ============================================================================

/// Exact analytical constant Omega_a for PR EOS:
/// Omega_a = (1 + (4 - 2*sqrt(2))^(1/3) + (4 + 2*sqrt(2))^(1/3)) / ... ~ 0.4572355289...
inline constexpr double OMEGA_A = 0.4572355289213822;

/// Exact analytical constant Omega_b for PR EOS:
/// Omega_b = (sqrt(2) - 1) / (3 * (1 + sqrt(2))) ~ 0.0777960739...
inline constexpr double OMEGA_B = 0.07779607390388846;

/// Square root of 2 mathematical constant
inline constexpr double SQRT2 = 1.4142135623730950488;

// ============================================================================
// Numerical Tolerances and Machine Precision Thresholds
// ============================================================================

/// Maximum acceptable absolute polynomial residual |f(Z)| for cubic roots
inline constexpr double CUBIC_POLY_RESIDUAL_TOL = 1.0e-12;

/// Tolerance below which cubic polynomial discriminant Delta is treated as zero
inline constexpr double CUBIC_DISCRIMINANT_ZERO_TOL = 1.0e-14;

/// Distance threshold to merge repeated/duplicate roots
inline constexpr double CUBIC_ROOT_MERGE_TOL = 1.0e-9;

/// Relative tolerance for PVT reciprocity: |P(T, v(T, P)) - P| / P (for P >= 1e-2 Pa)
inline constexpr double PVT_RECIPROCITY_TOL = 1.0e-10;

/// Relative tolerance for PVT reciprocity at ultra-low pressure (P < 1e-2 Pa)
inline constexpr double PVT_RECIPROCITY_LOW_P_TOL = 1.0e-6;

/// Relative tolerance for analytical vs numerical central-difference derivative da/dT
inline constexpr double DERIVATIVE_REL_TOL = 1.0e-7;

/// Dimensionless threshold for PR parameter B below which 3-term Taylor expansion is used
/// for (1 / (2*sqrt(2)*B)) * ln(...) to avoid catastrophic cancellation in floating-point 0/0
inline constexpr double SMALL_B_TAYLOR_THRESHOLD = 1.0e-5;

/// Relative tolerance on fugacity coefficients to identify liquid-vapor saturation equilibrium
inline constexpr double EQUILIBRIUM_FUGACITY_TOL = 1.0e-5;

} // namespace thermo::constants

#endif // THERMO_CONSTANTS_HPP
