#include "thermo/PengRobinsonEOS.hpp"
#include "thermo/Constants.hpp"
#include "thermo/CubicSolver.hpp"

#include <cmath>
#include <sstream>
#include <iomanip>
#include <stdexcept>
#include <algorithm>

namespace thermo {

void PengRobinsonEOS::validate_finite(double val, const std::string& name) {
    if (!std::isfinite(val)) {
        throw std::invalid_argument(name + " must be a finite numerical value.");
    }
}

PengRobinsonEOS::PengRobinsonEOS(PureComponent component)
    : component_(std::move(component)) {
    component_.validate();

    // Critical PR parameters
    a_c_ = constants::OMEGA_A * (constants::R * constants::R) * (component_.Tc * component_.Tc) / component_.Pc;
    b_ = constants::OMEGA_B * constants::R * component_.Tc / component_.Pc;
    kappa_ = 0.37464 + 1.54226 * component_.omega - 0.26992 * (component_.omega * component_.omega);

    init_reference_departures();
}

double PengRobinsonEOS::Tr(double T) const {
    validate_finite(T, "Temperature");
    if (T <= 0.0) {
        throw std::invalid_argument("Absolute temperature T must be strictly positive (in K).");
    }
    return T / component_.Tc;
}

double PengRobinsonEOS::alpha(double T) const {
    const double tr = Tr(T);
    const double sqrt_tr = std::sqrt(tr);
    const double term = 1.0 + kappa_ * (1.0 - sqrt_tr);
    return term * term;
}

double PengRobinsonEOS::a(double T) const {
    return a_c_ * alpha(T);
}

double PengRobinsonEOS::da_dT(double T) const {
    validate_finite(T, "Temperature");
    if (T <= 0.0) {
        throw std::invalid_argument("Absolute temperature T must be strictly positive (in K).");
    }

    const double tr = T / component_.Tc;
    const double sqrt_tr = std::sqrt(tr);
    const double term = 1.0 + kappa_ * (1.0 - sqrt_tr);

    // Exact analytical derivative:
    // da/dT = -a_c * kappa * [1 + kappa*(1 - sqrt(Tr))] / sqrt(T * Tc)
    const double denom = std::sqrt(T * component_.Tc);
    return -a_c_ * kappa_ * term / denom;
}

double PengRobinsonEOS::A_parameter(double T, double P) const {
    validate_finite(T, "Temperature");
    validate_finite(P, "Pressure");
    if (T <= 0.0) {
        throw std::invalid_argument("Absolute temperature T must be strictly positive (in K).");
    }
    if (P <= 0.0) {
        throw std::invalid_argument("Absolute pressure P must be strictly positive (in Pa).");
    }

    const double rt = constants::R * T;
    return a(T) * P / (rt * rt);
}

double PengRobinsonEOS::B_parameter(double T, double P) const {
    validate_finite(T, "Temperature");
    validate_finite(P, "Pressure");
    if (T <= 0.0) {
        throw std::invalid_argument("Absolute temperature T must be strictly positive (in K).");
    }
    if (P <= 0.0) {
        throw std::invalid_argument("Absolute pressure P must be strictly positive (in Pa).");
    }

    return b_ * P / (constants::R * T);
}

EOSParameters PengRobinsonEOS::parameters(double T) const {
    validate_finite(T, "Temperature");
    if (T <= 0.0) {
        throw std::invalid_argument("Absolute temperature T must be strictly positive (in K).");
    }

    EOSParameters p;
    p.Tr = Tr(T);
    p.alpha = alpha(T);
    p.a = a(T);
    p.da_dT = da_dT(T);
    p.ac = a_c_;
    p.b = b_;
    p.kappa = kappa_;
    return p;
}

double PengRobinsonEOS::calculate_pressure(double T, double v) const {
    validate_finite(T, "Temperature");
    validate_finite(v, "Molar volume");
    if (T <= 0.0) {
        throw std::invalid_argument("Absolute temperature T must be strictly positive (in K).");
    }
    if (v <= b_) {
        std::ostringstream oss;
        oss << "Molar volume v = " << std::scientific << std::setprecision(6) << v
            << " m^3/mol must be strictly greater than co-volume b = " << b_ << " m^3/mol.";
        throw std::invalid_argument(oss.str());
    }

    // PR attractive denominator: v^2 + 2*b*v - b^2
    const double denom_att = v * v + 2.0 * b_ * v - b_ * b_;
    if (denom_att <= 0.0) {
        throw std::invalid_argument("Molar volume violates Peng-Robinson attractive denominator domain.");
    }

    const double P = (constants::R * T) / (v - b_) - a(T) / denom_att;
    validate_finite(P, "Calculated pressure");
    return P;
}

double PengRobinsonEOS::calculate_dP_dv(double T, double v) const {
    validate_finite(T, "Temperature");
    validate_finite(v, "Molar volume");
    if (T <= 0.0) {
        throw std::invalid_argument("Absolute temperature T must be strictly positive (in K).");
    }
    if (v <= b_) {
        std::ostringstream oss;
        oss << "Molar volume v = " << std::scientific << std::setprecision(6) << v
            << " m^3/mol must be strictly greater than co-volume b = " << b_ << " m^3/mol.";
        throw std::invalid_argument(oss.str());
    }

    const double v_minus_b = v - b_;
    const double denom_att = v * v + 2.0 * b_ * v - b_ * b_;

    // Analytical isothermal derivative:
    // (dP/dv)_T = -R*T / (v - b)^2 + a(T)*(2v + 2b) / (v^2 + 2bv - b^2)^2
    const double term_rep = -(constants::R * T) / (v_minus_b * v_minus_b);
    const double term_att = a(T) * (2.0 * v + 2.0 * b_) / (denom_att * denom_att);

    return term_rep + term_att;
}

PVTResult PengRobinsonEOS::calculate_pvt(double T, double P) const {
    validate_finite(T, "Temperature");
    validate_finite(P, "Pressure");
    if (T <= 0.0) {
        throw std::invalid_argument("Absolute temperature T must be strictly positive (in K).");
    }
    if (P <= 0.0) {
        throw std::invalid_argument("Absolute pressure P must be strictly positive (in Pa).");
    }

    const double A = A_parameter(T, P);
    const double B = B_parameter(T, P);

    // Dimensionless PR cubic polynomial coefficients:
    // Z^3 - (1 - B)*Z^2 + (A - 3*B^2 - 2*B)*Z - (A*B - B^2 - B^3) = 0
    const double B2 = B * B;
    const double B3 = B2 * B;

    const double c2 = -(1.0 - B);
    const double c1 = A - 3.0 * B2 - 2.0 * B;
    const double c0 = -(A * B - B2 - B3);

    const CubicSolution cubic_sol = CubicSolver::solve(c2, c1, c0);

    PVTResult result;
    result.T = T;
    result.P = P;
    result.A = A;
    result.B = B;
    result.discriminant = cubic_sol.discriminant;

    const bool has_shomate = component_.shomate.has_value() &&
                             component_.shomate->is_defined() &&
                             component_.shomate->is_within_range(T) &&
                             component_.shomate->is_within_range(constants::T_REF);

    // Physical root filtering: retain only distinct roots with Z > B
    for (double Z : cubic_sol.distinct_roots) {
        if (Z > (B + 1.0e-12)) {
            const double v = Z * constants::R * T / P;
            const double dP_dv = calculate_dP_dv(T, v);

            PVTRoot root;
            root.Z = Z;
            root.molar_volume = v;
            root.poly_residual = std::abs(CubicSolver::evaluate_poly(Z, c2, c1, c0));
            root.dP_dv = dP_dv;
            root.stability = (dP_dv < 0.0) ? MechanicalStability::STABLE : MechanicalStability::UNSTABLE;
            root.is_physical = true;
            root.phase = PhaseTag::SINGLE_PHASE; // will be classified below

            // Session 4 Thermodynamic Properties:
            root.H_res = calculate_H_res(T, P, Z);
            root.S_res = calculate_S_res(T, P, Z);
            root.ln_phi = calculate_ln_phi(T, P, Z);
            root.phi = std::exp(root.ln_phi);
            root.fugacity = root.phi * P;

            if (has_shomate) {
                root.H_total = calculate_H_total(T, P, Z);
                root.S_total = calculate_S_total(T, P, Z);
            } else {
                root.H_total = std::nullopt;
                root.S_total = std::nullopt;
            }

            result.roots.push_back(root);
        }
    }

    // Sort physical roots in ascending order (Z_min <= ... <= Z_max)
    std::sort(result.roots.begin(), result.roots.end(),
              [](const PVTRoot& r1, const PVTRoot& r2) { return r1.Z < r2.Z; });

    // Classify phases according to root count and ordering
    if (result.roots.size() == 1) {
        result.roots[0].phase = PhaseTag::SINGLE_PHASE;
    } else if (result.roots.size() == 2) {
        result.roots[0].phase = PhaseTag::LIQUID;
        result.roots[1].phase = PhaseTag::VAPOR;
    } else if (result.roots.size() >= 3) {
        result.roots.front().phase = PhaseTag::LIQUID;
        result.roots.back().phase = PhaseTag::VAPOR;
        for (size_t i = 1; i + 1 < result.roots.size(); ++i) {
            result.roots[i].phase = PhaseTag::INTERMEDIATE;
        }
    }

    // Thermodynamic phase equilibrium and preference analysis
    result.equilibrium = evaluate_phase_preference(result);

    return result;
}

void PengRobinsonEOS::init_reference_departures() {
    const double A_ref = A_parameter(constants::T_REF, constants::P_REF);
    const double B_ref = B_parameter(constants::T_REF, constants::P_REF);

    const double B2 = B_ref * B_ref;
    const double B3 = B2 * B_ref;

    const double c2 = -(1.0 - B_ref);
    const double c1 = A_ref - 3.0 * B2 - 2.0 * B_ref;
    const double c0 = -(A_ref * B_ref - B2 - B3);

    const CubicSolution sol = CubicSolver::solve(c2, c1, c0);

    // Physical root filtering: find the vapor-like root (largest Z > B)
    double Z_vap_ref = 0.0;
    bool found = false;
    for (auto it = sol.distinct_roots.rbegin(); it != sol.distinct_roots.rend(); ++it) {
        if (*it > (B_ref + 1.0e-12)) {
            Z_vap_ref = *it;
            found = true;
            break;
        }
    }

    if (found) {
        H_res_ref_ = calculate_H_res(constants::T_REF, constants::P_REF, Z_vap_ref);
        S_res_ref_ = calculate_S_res(constants::T_REF, constants::P_REF, Z_vap_ref);
        ref_state_initialized_ = true;
    }
}

std::pair<double, double> PengRobinsonEOS::reference_state_departures() const {
    if (!ref_state_initialized_) {
        throw std::runtime_error("Reference state departures have not been initialized.");
    }
    return {H_res_ref_, S_res_ref_};
}

double PengRobinsonEOS::calculate_I_L(double Z, double B) const {
    validate_finite(Z, "Compressibility factor Z");
    validate_finite(B, "Dimensionless parameter B");

    if (Z <= B) {
        std::ostringstream oss;
        oss << "Compressibility factor Z = " << std::fixed << std::setprecision(6) << Z
            << " must be strictly greater than parameter B = " << B << ".";
        throw std::invalid_argument(oss.str());
    }

    const double denom = Z + (1.0 - constants::SQRT2) * B;
    if (denom <= 0.0) {
        throw std::invalid_argument("Compressibility factor violates PR logarithmic denominator domain.");
    }

    const double u = 2.0 * constants::SQRT2 * B / denom;

    if (std::abs(u) < constants::SMALL_B_TAYLOR_THRESHOLD) {
        // 4-term Taylor expansion for ln(1+u)/u:
        // ln(1+u)/u = 1 - u/2 + u^2/3 - u^3/4
        const double series = 1.0 - 0.5 * u + (1.0 / 3.0) * u * u - 0.25 * u * u * u;
        return (1.0 / denom) * series;
    }

    const double L = std::log1p(u);
    return L / (2.0 * constants::SQRT2 * B);
}

double PengRobinsonEOS::calculate_L(double Z, double B) const {
    const double I_L = calculate_I_L(Z, B);
    return 2.0 * constants::SQRT2 * B * I_L;
}

double PengRobinsonEOS::calculate_H_res(double T, double P, double Z) const {
    validate_finite(T, "Temperature");
    validate_finite(P, "Pressure");
    validate_finite(Z, "Compressibility factor Z");
    if (T <= 0.0) {
        throw std::invalid_argument("Absolute temperature T must be strictly positive (in K).");
    }
    if (P <= 0.0) {
        throw std::invalid_argument("Absolute pressure P must be strictly positive (in Pa).");
    }

    const double B = B_parameter(T, P);
    if (Z <= B) {
        std::ostringstream oss;
        oss << "Compressibility factor Z = " << std::fixed << std::setprecision(6) << Z
            << " must be strictly greater than parameter B = " << B << ".";
        throw std::invalid_argument(oss.str());
    }

    const double L = calculate_L(Z, B);
    const double da = da_dT(T);
    const double a_val = a(T);

    // H^res = R*T*(Z - 1) + ((T*da/dT - a) / (2*sqrt(2)*b)) * L
    const double term1 = constants::R * T * (Z - 1.0);
    const double term2 = ((T * da - a_val) / (2.0 * constants::SQRT2 * b_)) * L;

    const double H_res = term1 + term2;
    validate_finite(H_res, "Residual enthalpy");
    return H_res;
}

double PengRobinsonEOS::calculate_S_res(double T, double P, double Z) const {
    validate_finite(T, "Temperature");
    validate_finite(P, "Pressure");
    validate_finite(Z, "Compressibility factor Z");
    if (T <= 0.0) {
        throw std::invalid_argument("Absolute temperature T must be strictly positive (in K).");
    }
    if (P <= 0.0) {
        throw std::invalid_argument("Absolute pressure P must be strictly positive (in Pa).");
    }

    const double B = B_parameter(T, P);
    if (Z <= B) {
        std::ostringstream oss;
        oss << "Compressibility factor Z = " << std::fixed << std::setprecision(6) << Z
            << " must be strictly greater than parameter B = " << B << ".";
        throw std::invalid_argument(oss.str());
    }

    const double L = calculate_L(Z, B);
    const double da = da_dT(T);

    // S^res = R*ln(Z - B) + (da/dT / (2*sqrt(2)*b)) * L
    const double eps = Z - 1.0 - B;
    const double log_Z_minus_B = (std::abs(eps) < constants::SMALL_B_TAYLOR_THRESHOLD)
                                     ? std::log1p(eps)
                                     : std::log(Z - B);

    const double term1 = constants::R * log_Z_minus_B;
    const double term2 = (da / (2.0 * constants::SQRT2 * b_)) * L;

    const double S_res = term1 + term2;
    validate_finite(S_res, "Residual entropy");
    return S_res;
}

double PengRobinsonEOS::calculate_ln_phi(double T, double P, double Z) const {
    validate_finite(T, "Temperature");
    validate_finite(P, "Pressure");
    validate_finite(Z, "Compressibility factor Z");
    if (T <= 0.0) {
        throw std::invalid_argument("Absolute temperature T must be strictly positive (in K).");
    }
    if (P <= 0.0) {
        throw std::invalid_argument("Absolute pressure P must be strictly positive (in Pa).");
    }

    const double B = B_parameter(T, P);
    if (Z <= B) {
        std::ostringstream oss;
        oss << "Compressibility factor Z = " << std::fixed << std::setprecision(6) << Z
            << " must be strictly greater than parameter B = " << B << ".";
        throw std::invalid_argument(oss.str());
    }

    const double A = A_parameter(T, P);
    const double I_L = calculate_I_L(Z, B);

    // Regularization of (Z - 1) - ln(Z - B)
    // eps = Z - 1 - B. ln(Z - B) = ln(1 + eps)
    // (Z - 1) - ln(Z - B) = B + eps^2 * (1/2 - eps/3 + eps^2/4)
    const double eps = Z - 1.0 - B;
    double log_diff = 0.0;
    if (std::abs(eps) < constants::SMALL_B_TAYLOR_THRESHOLD) {
        const double series = 0.5 - (1.0 / 3.0) * eps + 0.25 * eps * eps;
        log_diff = B + eps * eps * series;
    } else {
        log_diff = (Z - 1.0) - std::log(Z - B);
    }

    const double ln_phi = log_diff - A * I_L;
    validate_finite(ln_phi, "Log fugacity coefficient");
    return ln_phi;
}

double PengRobinsonEOS::calculate_phi(double T, double P, double Z) const {
    const double ln_phi = calculate_ln_phi(T, P, Z);
    const double phi = std::exp(ln_phi);
    validate_finite(phi, "Fugacity coefficient");
    return phi;
}

double PengRobinsonEOS::calculate_fugacity(double T, double P, double Z) const {
    const double phi = calculate_phi(T, P, Z);
    const double f = phi * P;
    validate_finite(f, "Fugacity");
    return f;
}

double PengRobinsonEOS::calculate_H_total(double T, double P, double Z) const {
    if (!component_.shomate.has_value() || !component_.shomate->is_defined()) {
        throw std::runtime_error("Component '" + component_.name +
                                 "' does not have ideal-gas Shomate heat capacity parameters defined.");
    }
    if (!component_.shomate->is_within_range(T)) {
        std::ostringstream oss;
        oss << "Temperature T = " << std::fixed << std::setprecision(2) << T
            << " K is outside valid Shomate interval [" << component_.shomate->T_min << ", "
            << component_.shomate->T_max << "] K for component '" << component_.name << "'.";
        throw std::out_of_range(oss.str());
    }
    if (!component_.shomate->is_within_range(constants::T_REF)) {
        std::ostringstream oss;
        oss << "Reference temperature T_ref = " << std::fixed << std::setprecision(2) << constants::T_REF
            << " K is outside valid Shomate interval [" << component_.shomate->T_min << ", "
            << component_.shomate->T_max << "] K for component '" << component_.name << "'.";
        throw std::out_of_range(oss.str());
    }

    const double dH_ig = component_.shomate->integrate_cp(constants::T_REF, T);
    const double H_res = calculate_H_res(T, P, Z);
    return dH_ig + H_res - H_res_ref_;
}

double PengRobinsonEOS::calculate_S_total(double T, double P, double Z) const {
    if (!component_.shomate.has_value() || !component_.shomate->is_defined()) {
        throw std::runtime_error("Component '" + component_.name +
                                 "' does not have ideal-gas Shomate heat capacity parameters defined.");
    }
    if (!component_.shomate->is_within_range(T)) {
        std::ostringstream oss;
        oss << "Temperature T = " << std::fixed << std::setprecision(2) << T
            << " K is outside valid Shomate interval [" << component_.shomate->T_min << ", "
            << component_.shomate->T_max << "] K for component '" << component_.name << "'.";
        throw std::out_of_range(oss.str());
    }
    if (!component_.shomate->is_within_range(constants::T_REF)) {
        std::ostringstream oss;
        oss << "Reference temperature T_ref = " << std::fixed << std::setprecision(2) << constants::T_REF
            << " K is outside valid Shomate interval [" << component_.shomate->T_min << ", "
            << component_.shomate->T_max << "] K for component '" << component_.name << "'.";
        throw std::out_of_range(oss.str());
    }

    const double dS_ig = component_.shomate->integrate_cp_over_T(constants::T_REF, T);
    const double S_res = calculate_S_res(T, P, Z);
    return dS_ig - constants::R * std::log(P / constants::P_REF) + S_res - S_res_ref_;
}

PhaseEquilibriumResult PengRobinsonEOS::evaluate_phase_preference(const PVTResult& pvt) const {
    PhaseEquilibriumResult result;
    if (pvt.roots.empty()) {
        result.preference = PhasePreference::SINGLE_PHASE;
        result.summary = "No physical roots found.";
        return result;
    }
    if (pvt.roots.size() == 1) {
        result.preference = PhasePreference::SINGLE_PHASE;
        result.phi_liq = pvt.roots[0].phi;
        result.phi_vap = pvt.roots[0].phi;
        result.phi_ratio = 1.0;
        result.delta_G_res = 0.0;
        result.is_equilibrium = false;
        result.summary = "Single physical root; no phase coexistence.";
        return result;
    }

    const PVTRoot* liq = pvt.liquid_root();
    const PVTRoot* vap = pvt.vapor_root();

    if (!liq || !vap || liq == vap) {
        result.preference = PhasePreference::SINGLE_PHASE;
        result.summary = "No distinct liquid and vapor roots found.";
        return result;
    }

    result.phi_liq = liq->phi;
    result.phi_vap = vap->phi;
    result.phi_ratio = liq->phi / vap->phi;

    // Delta G_res = G_res,liq - G_res,vap = RT * ln(phi_liq / phi_vap)
    const double ln_ratio = std::log(result.phi_ratio);
    result.delta_G_res = constants::R * pvt.T * ln_ratio;

    // Relative difference in fugacity coefficients
    const double rel_diff = std::abs(liq->phi - vap->phi) / vap->phi;

    if (rel_diff <= constants::EQUILIBRIUM_FUGACITY_TOL || std::abs(ln_ratio) <= constants::EQUILIBRIUM_FUGACITY_TOL) {
        result.preference = PhasePreference::EQUILIBRIUM;
        result.is_equilibrium = true;
        result.summary = "Liquid and vapor phases are at equilibrium (saturation state).";
    } else if (liq->phi < vap->phi) {
        result.preference = PhasePreference::LIQUID_PREFERRED;
        result.is_equilibrium = false;
        result.summary = "Liquid phase preferred (lower fugacity / chemical potential; subcooled liquid).";
    } else {
        result.preference = PhasePreference::VAPOR_PREFERRED;
        result.is_equilibrium = false;
        result.summary = "Vapor phase preferred (lower fugacity / chemical potential; superheated vapor).";
    }

    return result;
}

} // namespace thermo
