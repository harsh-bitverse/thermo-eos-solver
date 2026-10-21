#ifndef THERMO_PENG_ROBINSON_EOS_HPP
#define THERMO_PENG_ROBINSON_EOS_HPP

#include "thermo/Component.hpp"
#include "thermo/StateProperties.hpp"
#include <string>

namespace thermo {

/**
 * @brief Peng-Robinson (1976) Equation of State implementation for pure components.
 * 
 * Pressure equation:
 *   P = R*T / (v - b) - a(T) / (v^2 + 2*b*v - b^2)
 * 
 * Dimensionless compressibility polynomial:
 *   Z^3 - (1 - B)*Z^2 + (A - 3*B^2 - 2*B)*Z - (A*B - B^2 - B^3) = 0
 * 
 * Features:
 * - Domain validation: strict rejection of T <= 0, P <= 0, v <= b, and non-finite numbers.
 * - Accurate precomputation of PR-1976 parameters (a_c, b, kappa).
 * - Exact analytical temperature derivative da/dT.
 * - Forward PVT: solves cubic in Z, discards unphysical roots (Z <= B), calculates v,
 *   evaluates mechanical stability (dP/dv)_T, and tags roots as VAPOR, LIQUID, or INTERMEDIATE.
 * - Inverse PVT: computes P directly from (T, v).
 * - Evaluation of isothermal derivative (dP/dv)_T.
 */
class PengRobinsonEOS {
public:
    /**
     * @brief Constructs a PengRobinsonEOS model for the given pure component.
     * @param component Validated pure component.
     * @throws std::invalid_argument if component parameters are non-physical.
     */
    explicit PengRobinsonEOS(PureComponent component);

    /// Returns the underlying pure component
    [[nodiscard]] const PureComponent& component() const noexcept { return component_; }

    /// Critical attractive parameter a_c in [Pa * m^6 / mol^2]
    [[nodiscard]] double ac() const noexcept { return a_c_; }

    /// Co-volume parameter b in [m^3 / mol]
    [[nodiscard]] double b() const noexcept { return b_; }

    /// Acentric factor parameter kappa
    [[nodiscard]] double kappa() const noexcept { return kappa_; }

    /// Reduced temperature T_r = T / T_c
    [[nodiscard]] double Tr(double T) const;

    /// PR alpha function alpha(T_r, omega)
    [[nodiscard]] double alpha(double T) const;

    /// Temperature-dependent parameter a(T) = a_c * alpha(T) in [Pa * m^6 / mol^2]
    [[nodiscard]] double a(double T) const;

    /// Exact analytical temperature derivative da/dT in [Pa * m^6 / (mol^2 * K)]
    [[nodiscard]] double da_dT(double T) const;

    /// Dimensionless attractive parameter A = a(T)*P / (R*T)^2
    [[nodiscard]] double A_parameter(double T, double P) const;

    /// Dimensionless repulsive parameter B = b*P / (R*T)
    [[nodiscard]] double B_parameter(double T, double P) const;

    /// Returns all EOS parameters evaluated at temperature T
    [[nodiscard]] EOSParameters parameters(double T) const;

    /**
     * @brief Forward PVT calculation: solves for compressibility factor Z and molar volume v at (T, P).
     * Retains only physical roots (Z > B), calculates (dP/dv)_T, and classifies roots.
     * @param T Absolute temperature in [K] (must be strictly positive).
     * @param P Absolute pressure in [Pa] (must be strictly positive).
     * @return PVTResult containing all physically admissible roots.
     * @throws std::invalid_argument if T <= 0, P <= 0, or non-finite.
     */
    [[nodiscard]] PVTResult calculate_pvt(double T, double P) const;

    /**
     * @brief Inverse PVT calculation: computes pressure P given temperature T and molar volume v.
     * Evaluates directly from the Peng-Robinson equation of state.
     * @param T Absolute temperature in [K] (must be strictly positive).
     * @param v Molar volume in [m^3 / mol] (must be strictly greater than b).
     * @return Calculated pressure in [Pa].
     * @throws std::invalid_argument if T <= 0, v <= b, or non-finite.
     */
    [[nodiscard]] double calculate_pressure(double T, double v) const;

    /**
     * @brief Evaluates the exact analytical isothermal derivative (dP/dv)_T in [Pa / (m^3 / mol)].
     * Formula: (dP/dv)_T = -R*T / (v - b)^2 + a(T)*(2v + 2b) / (v^2 + 2bv - b^2)^2
     * @param T Absolute temperature in [K] (must be strictly positive).
     * @param v Molar volume in [m^3 / mol] (must be strictly greater than b).
     * @throws std::invalid_argument if T <= 0, v <= b, or non-finite.
     */
    [[nodiscard]] double calculate_dP_dv(double T, double v) const;

private:
    PureComponent component_;
    double a_c_{0.0};
    double b_{0.0};
    double kappa_{0.0};

    static void validate_finite(double val, const std::string& name);
};

} // namespace thermo

#endif // THERMO_PENG_ROBINSON_EOS_HPP
