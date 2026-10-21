#ifndef THERMO_STATE_PROPERTIES_HPP
#define THERMO_STATE_PROPERTIES_HPP

#include <vector>
#include <string>
#include <optional>

namespace thermo {

/**
 * @brief Qualitative classification of compressibility-factor roots.
 */
enum class PhaseTag {
    SINGLE_PHASE, ///< Exactly one physical root exists (e.g. supercritical fluid or outside spinodal)
    VAPOR,        ///< Vapor-like root (largest physical compressibility factor)
    LIQUID,       ///< Liquid-like root (smallest physical compressibility factor)
    INTERMEDIATE  ///< Unphysical/metastable intermediate root in three-root region
};

/**
 * @brief Mechanical stability indicator based on the isothermal derivative (dP/dv)_T.
 */
enum class MechanicalStability {
    STABLE,  ///< (dP/dv)_T < 0: Mechanically stable (compression increases pressure)
    UNSTABLE ///< (dP/dv)_T >= 0: Mechanically unstable (spinodal unobservable state)
};

/**
 * @brief Thermodynamic phase preference based on chemical potential / fugacity comparison.
 */
enum class PhasePreference {
    SINGLE_PHASE,     ///< Exactly one physical root exists; no phase coexistence
    LIQUID_PREFERRED, ///< Multiple roots: phi_liq < phi_vap (G_liq < G_vap; subcooled liquid)
    VAPOR_PREFERRED,  ///< Multiple roots: phi_vap < phi_liq (G_vap < G_liq; superheated vapor)
    EQUILIBRIUM       ///< Multiple roots: phi_liq == phi_vap within tolerance (saturation coexistence)
};

/**
 * @brief Summary of pure-component phase equilibrium analysis for multi-root states.
 */
struct PhaseEquilibriumResult {
    PhasePreference preference{PhasePreference::SINGLE_PHASE};
    double delta_G_res{0.0};  ///< Residual Gibbs difference: G_res,liq - G_res,vap = RT*ln(phi_liq / phi_vap) in [J / mol]
    double phi_liq{0.0};      ///< Fugacity coefficient of liquid-like root
    double phi_vap{0.0};      ///< Fugacity coefficient of vapor-like root
    double phi_ratio{1.0};    ///< Fugacity coefficient ratio phi_liq / phi_vap
    bool is_equilibrium{false};///< True if |phi_liq - phi_vap| / phi_vap <= EQUILIBRIUM_FUGACITY_TOL
    std::string summary;      ///< Human-readable diagnostic description
};

/**
 * @brief Representation of a single compressibility-factor root and its thermodynamic state properties.
 */
struct PVTRoot {
    double Z{1.0};               ///< Compressibility factor Z = P*v / (R*T)
    double molar_volume{0.0};    ///< Molar volume v in [m^3 / mol]
    double poly_residual{0.0};   ///< Absolute residual |P(Z)| of the cubic polynomial
    double dP_dv{0.0};           ///< Isothermal derivative (dP/dv)_T in [Pa / (m^3 / mol)]
    MechanicalStability stability{MechanicalStability::STABLE};
    PhaseTag phase{PhaseTag::SINGLE_PHASE};
    bool is_physical{false};     ///< True if Z > B (molar volume v > co-volume b)

    // Thermodynamic properties for this root:
    double H_res{0.0};           ///< Residual enthalpy in [J / mol]
    double S_res{0.0};           ///< Residual entropy in [J / (mol * K)]
    double ln_phi{0.0};          ///< Natural log of fugacity coefficient [dimensionless]
    double phi{1.0};             ///< Fugacity coefficient [dimensionless]
    double fugacity{0.0};        ///< Fugacity f = phi * P in [Pa]
    std::optional<double> H_total{std::nullopt}; ///< Reference-state total enthalpy in [J / mol]
    std::optional<double> S_total{std::nullopt}; ///< Reference-state total entropy in [J / (mol * K)]
};

/**
 * @brief Complete result container for forward PVT and thermodynamic calculations at a given (T, P).
 */
struct PVTResult {
    double T{0.0};               ///< Temperature in [K]
    double P{0.0};               ///< Pressure in [Pa]
    double A{0.0};               ///< Dimensionless PR attractive parameter a*P / (R*T)^2
    double B{0.0};               ///< Dimensionless PR repulsive parameter b*P / (R*T)
    double discriminant{0.0};    ///< Depressed cubic discriminant Delta
    std::vector<PVTRoot> roots;  ///< Retained physically admissible roots (Z > B), sorted ascending
    PhaseEquilibriumResult equilibrium; ///< Thermodynamic phase equilibrium analysis

    /// Returns true if multiple physically admissible roots exist at this state
    [[nodiscard]] bool has_multiple_roots() const noexcept {
        return roots.size() > 1;
    }

    /// Returns pointer to vapor-like root (largest Z), or nullptr if no roots exist
    [[nodiscard]] const PVTRoot* vapor_root() const noexcept {
        for (auto it = roots.rbegin(); it != roots.rend(); ++it) {
            if (it->phase == PhaseTag::VAPOR || it->phase == PhaseTag::SINGLE_PHASE) {
                return &(*it);
            }
        }
        return roots.empty() ? nullptr : &roots.back();
    }

    /// Returns pointer to liquid-like root (smallest Z), or nullptr if not present
    [[nodiscard]] const PVTRoot* liquid_root() const noexcept {
        for (const auto& r : roots) {
            if (r.phase == PhaseTag::LIQUID) {
                return &r;
            }
        }
        return nullptr;
    }

    /// Returns pointer to intermediate root (unstable), or nullptr if not present
    [[nodiscard]] const PVTRoot* intermediate_root() const noexcept {
        for (const auto& r : roots) {
            if (r.phase == PhaseTag::INTERMEDIATE) {
                return &r;
            }
        }
        return nullptr;
    }

    /// Returns pointer to the thermodynamically preferred root based on equilibrium analysis
    [[nodiscard]] const PVTRoot* preferred_root() const noexcept {
        if (roots.empty()) return nullptr;
        if (roots.size() == 1) return &roots[0];
        if (equilibrium.preference == PhasePreference::LIQUID_PREFERRED) {
            return liquid_root();
        }
        return vapor_root();
    }
};

/**
 * @brief Complete temperature-dependent parameter evaluation for Peng-Robinson EOS.
 */
struct EOSParameters {
    double a{0.0};     ///< Temperature-dependent attractive parameter a(T) [Pa * m^6 / mol^2]
    double b{0.0};     ///< Co-volume parameter b [m^3 / mol]
    double da_dT{0.0}; ///< Analytical derivative da/dT [Pa * m^6 / (mol^2 * K)]
    double ac{0.0};    ///< Critical attractive parameter a_c [Pa * m^6 / mol^2]
    double kappa{0.0}; ///< PR acentric factor parameter
    double alpha{0.0}; ///< PR alpha function [dimensionless]
    double Tr{0.0};    ///< Reduced temperature T / T_c
};

} // namespace thermo

#endif // THERMO_STATE_PROPERTIES_HPP
