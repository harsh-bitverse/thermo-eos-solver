#ifndef THERMO_PENG_ROBINSON_EOS_HPP
#define THERMO_PENG_ROBINSON_EOS_HPP

#include "thermo/Component.hpp"
#include "thermo/StateProperties.hpp"
#include <string>
#include <utility>

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

    // ========================================================================
    // Session 4: Residual Properties, Fugacity, Total H/S, & Phase Equilibrium
    // ========================================================================

    /**
     * @brief Dimensionless regularized logarithmic ratio: I_L = L / (2*sqrt(2)*B).
     * Controlled via 4-term Taylor expansion below SMALL_B_TAYLOR_THRESHOLD to prevent 0/0 breakdown.
     * @param Z Compressibility factor (must satisfy Z > B).
     * @param B Dimensionless PR parameter b*P / (R*T).
     */
    [[nodiscard]] double calculate_I_L(double Z, double B) const;

    /**
     * @brief PR logarithmic departure ratio: L = ln((Z + (1+sqrt2)B) / (Z + (1-sqrt2)B)).
     * Evaluated using std::log1p and regularized below SMALL_B_TAYLOR_THRESHOLD.
     */
    [[nodiscard]] double calculate_L(double Z, double B) const;

    /**
     * @brief Evaluates residual enthalpy H^res(T, P, Z) in [J / mol].
     * Formula: H^res = R*T*(Z - 1) + ((T*da/dT - a) / (2*sqrt(2)*b)) * L
     * @throws std::invalid_argument if inputs are invalid or Z <= B.
     */
    [[nodiscard]] double calculate_H_res(double T, double P, double Z) const;

    /**
     * @brief Evaluates residual entropy S^res(T, P, Z) in [J / (mol * K)].
     * Formula: S^res = R*ln(Z - B) + (da/dT / (2*sqrt(2)*b)) * L
     * @throws std::invalid_argument if inputs are invalid or Z <= B.
     */
    [[nodiscard]] double calculate_S_res(double T, double P, double Z) const;

    /**
     * @brief Evaluates natural logarithm of pure-component fugacity coefficient ln(phi).
     * Formula: ln(phi) = (Z - 1) - ln(Z - B) - (A / (2*sqrt(2)*B)) * L
     * Uses small-B / small-epsilon Taylor series below SMALL_B_TAYLOR_THRESHOLD.
     * @throws std::invalid_argument if inputs are invalid or Z <= B.
     */
    [[nodiscard]] double calculate_ln_phi(double T, double P, double Z) const;

    /**
     * @brief Evaluates pure-component fugacity coefficient phi = exp(ln(phi)) [dimensionless].
     */
    [[nodiscard]] double calculate_phi(double T, double P, double Z) const;

    /**
     * @brief Evaluates pure-component fugacity f = phi * P in [Pa].
     */
    [[nodiscard]] double calculate_fugacity(double T, double P, double Z) const;

    /**
     * @brief Evaluates reference-state total enthalpy H(T, P, Z) in [J / mol] relative to
     * T_ref = 298.15 K and P_ref = 100,000 Pa.
     * Formula: H = Delta H^ig(T_ref -> T) + H^res(T, P, Z) - H^res(T_ref, P_ref, Z_ref)
     * @throws std::runtime_error if component does not have Shomate parameters.
     * @throws std::out_of_range if T or T_ref is outside the valid Shomate interval.
     */
    [[nodiscard]] double calculate_H_total(double T, double P, double Z) const;

    /**
     * @brief Evaluates reference-state total entropy S(T, P, Z) in [J / (mol * K)] relative to
     * T_ref = 298.15 K and P_ref = 100,000 Pa.
     * Formula: S = Delta S^ig(T_ref -> T) - R*ln(P / P_ref) + S^res(T, P, Z) - S^res(T_ref, P_ref, Z_ref)
     * @throws std::runtime_error if component does not have Shomate parameters.
     * @throws std::out_of_range if T or T_ref is outside the valid Shomate interval.
     */
    [[nodiscard]] double calculate_S_total(double T, double P, double Z) const;

    /**
     * @brief Evaluates thermodynamic phase preference for a multi-root PVTResult.
     * Compares liquid vs vapor fugacities and excludes mechanically unstable intermediate roots.
     */
    [[nodiscard]] PhaseEquilibriumResult evaluate_phase_preference(const PVTResult& pvt) const;

    /**
     * @brief Returns reference-state departures (H_res_ref, S_res_ref) at 298.15 K, 100,000 Pa.
     */
    [[nodiscard]] std::pair<double, double> reference_state_departures() const;

private:
    PureComponent component_;
    double a_c_{0.0};
    double b_{0.0};
    double kappa_{0.0};

    // Reference state departure cache (T_ref = 298.15 K, P_ref = 100,000 Pa, vapor-like root)
    double H_res_ref_{0.0};
    double S_res_ref_{0.0};
    bool ref_state_initialized_{false};

    void init_reference_departures();
    static void validate_finite(double val, const std::string& name);
};

} // namespace thermo

#endif // THERMO_PENG_ROBINSON_EOS_HPP
