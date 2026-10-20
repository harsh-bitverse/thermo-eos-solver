#include "thermo/Component.hpp"

#include <cmath>
#include <sstream>
#include <iomanip>

namespace thermo {

bool ShomateParameters::is_within_range(double T) const noexcept {
    // Add small epsilon tolerance of 1e-4 K to handle floating-point edge cases
    return (T >= (T_min - 1.0e-4)) && (T <= (T_max + 1.0e-4));
}

bool ShomateParameters::is_defined() const noexcept {
    return (T_min > 0.0) && (T_max > T_min);
}

double ShomateParameters::evaluate_cp(double T) const {
    if (!is_within_range(T)) {
        std::ostringstream oss;
        oss << "Temperature T = " << std::fixed << std::setprecision(2) << T
            << " K is outside valid NIST Shomate interval [" << T_min << ", " << T_max << "] K.";
        throw std::out_of_range(oss.str());
    }

    const double t = T / 1000.0;
    const double t2 = t * t;
    const double t3 = t2 * t;

    return A + B * t + C * t2 + D * t3 + (E / t2);
}

double ShomateParameters::integrate_cp(double T_ref, double T) const {
    if (!is_within_range(T)) {
        std::ostringstream oss;
        oss << "Target temperature T = " << std::fixed << std::setprecision(2) << T
            << " K is outside valid NIST Shomate interval [" << T_min << ", " << T_max << "] K.";
        throw std::out_of_range(oss.str());
    }
    if (!is_within_range(T_ref)) {
        std::ostringstream oss;
        oss << "Reference temperature T_ref = " << std::fixed << std::setprecision(2) << T_ref
            << " K is outside valid NIST Shomate interval [" << T_min << ", " << T_max << "] K.";
        throw std::out_of_range(oss.str());
    }

    const double t = T / 1000.0;
    const double t0 = T_ref / 1000.0;

    const double t2 = t * t;
    const double t0_2 = t0 * t0;

    const double t3 = t2 * t;
    const double t0_3 = t0_2 * t0;

    const double t4 = t3 * t;
    const double t0_4 = t0_3 * t0;

    // Definite integral \int_{T_ref}^{T} Cp(T') dT'
    // = 1000 * [ A*(t - t0) + (B/2)*(t^2 - t0^2) + (C/3)*(t^3 - t0^3) + (D/4)*(t^4 - t0^4) - E*(1/t - 1/t0) ]
    const double term = A * (t - t0)
                      + 0.5 * B * (t2 - t0_2)
                      + (1.0 / 3.0) * C * (t3 - t0_3)
                      + 0.25 * D * (t4 - t0_4)
                      - E * ((1.0 / t) - (1.0 / t0));

    return 1000.0 * term; // Convert to J / mol
}

double ShomateParameters::integrate_cp_over_T(double T_ref, double T) const {
    if (!is_within_range(T)) {
        std::ostringstream oss;
        oss << "Target temperature T = " << std::fixed << std::setprecision(2) << T
            << " K is outside valid NIST Shomate interval [" << T_min << ", " << T_max << "] K.";
        throw std::out_of_range(oss.str());
    }
    if (!is_within_range(T_ref)) {
        std::ostringstream oss;
        oss << "Reference temperature T_ref = " << std::fixed << std::setprecision(2) << T_ref
            << " K is outside valid NIST Shomate interval [" << T_min << ", " << T_max << "] K.";
        throw std::out_of_range(oss.str());
    }

    const double t = T / 1000.0;
    const double t0 = T_ref / 1000.0;

    const double t2 = t * t;
    const double t0_2 = t0 * t0;

    const double t3 = t2 * t;
    const double t0_3 = t0_2 * t0;

    // Definite integral \int_{T_ref}^{T} (Cp(T') / T') dT'
    // Note: dT' / T' = dt' / t'.
    // Integral of A/t + B + C*t + D*t^2 + E/t^3
    // = A*ln(t/t0) + B*(t - t0) + (C/2)*(t^2 - t0^2) + (D/3)*(t^3 - t0^3) - (E/2)*(1/t^2 - 1/t0^2)
    const double term = A * std::log(t / t0)
                      + B * (t - t0)
                      + 0.5 * C * (t2 - t0_2)
                      + (1.0 / 3.0) * D * (t3 - t0_3)
                      - 0.5 * E * ((1.0 / t2) - (1.0 / t0_2));

    return term; // in J / (mol * K)
}

bool PureComponent::is_valid() const noexcept {
    return (!name.empty()) &&
           (molar_mass > 0.0) &&
           (Tc > 0.0) &&
           (Pc > 0.0) &&
           (omega > -1.0);
}

void PureComponent::validate() const {
    if (name.empty()) {
        throw std::invalid_argument("Component name cannot be empty.");
    }
    if (molar_mass <= 0.0) {
        throw std::invalid_argument("Molar mass must be strictly positive (in kg/mol).");
    }
    if (Tc <= 0.0) {
        throw std::invalid_argument("Critical temperature Tc must be strictly positive (in K).");
    }
    if (Pc <= 0.0) {
        throw std::invalid_argument("Critical pressure Pc must be strictly positive (in Pa).");
    }
    if (omega <= -1.0) {
        throw std::invalid_argument("Acentric factor omega must be greater than -1.0.");
    }
}

} // namespace thermo
