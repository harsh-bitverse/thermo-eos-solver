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

    return result;
}

} // namespace thermo
