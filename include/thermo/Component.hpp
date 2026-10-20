#ifndef THERMO_COMPONENT_HPP
#define THERMO_COMPONENT_HPP

#include <string>
#include <optional>
#include <stdexcept>

namespace thermo {

/**
 * @brief Parameters and analytical integrals for the authentic NIST Shomate equation.
 * 
 * Functional form (NIST Chemistry WebBook, SRD 69):
 *   t = T / 1000.0   (where T is in Kelvin)
 *   Cp(t) = A + B*t + C*t^2 + D*t^3 + E/(t^2)   [J / (mol * K)]
 * 
 * Note the inverse-square term E / t^2, which provides high fidelity across broad
 * temperature intervals and must not be confused with a standard polynomial.
 */
struct ShomateParameters {
    double A{0.0};
    double B{0.0};
    double C{0.0};
    double D{0.0};
    double E{0.0};
    double T_min{0.0}; // Valid lower temperature bound in [K]
    double T_max{0.0}; // Valid upper temperature bound in [K]

    /// Checks if a temperature T is within the validated range [T_min, T_max]
    [[nodiscard]] bool is_within_range(double T) const noexcept;

    /// Checks if this Shomate parameter set has been initialized with valid bounds
    [[nodiscard]] bool is_defined() const noexcept;

    /**
     * @brief Evaluates ideal-gas heat capacity Cp(T) in J / (mol * K).
     * @throws std::out_of_range if T is outside [T_min, T_max].
     */
    [[nodiscard]] double evaluate_cp(double T) const;

    /**
     * @brief Evaluates the definite integral \int_{T_ref}^{T} Cp(T') dT' in J / mol.
     * Evaluates to identically 0 when T == T_ref.
     * @throws std::out_of_range if T or T_ref is outside [T_min, T_max].
     */
    [[nodiscard]] double integrate_cp(double T_ref, double T) const;

    /**
     * @brief Evaluates the definite integral \int_{T_ref}^{T} (Cp(T') / T') dT' in J / (mol * K).
     * Evaluates to identically 0 when T == T_ref.
     * @throws std::out_of_range if T or T_ref is outside [T_min, T_max].
     */
    [[nodiscard]] double integrate_cp_over_T(double T_ref, double T) const;
};

/**
 * @brief Representation of a pure chemical component with critical properties and ideal-gas data.
 * All properties are stored strictly in SI units.
 */
struct PureComponent {
    std::string name;
    std::string formula;
    double molar_mass{0.0};  ///< Molar mass in [kg/mol] (e.g. 0.016043 for CH4)
    double Tc{0.0};          ///< Critical temperature in [K]
    double Pc{0.0};          ///< Critical pressure in [Pa]
    double omega{0.0};       ///< Pitzer acentric factor [dimensionless]
    std::optional<ShomateParameters> shomate{std::nullopt};
    std::string citation;    ///< Traceable data provenance

    /// Checks if component properties are physically valid (Tc > 0, Pc > 0, M > 0, omega > -1)
    [[nodiscard]] bool is_valid() const noexcept;

    /// Throws std::invalid_argument with a descriptive message if the component is invalid
    void validate() const;
};

} // namespace thermo

#endif // THERMO_COMPONENT_HPP
